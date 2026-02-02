import multiprocessing as _mp
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import polars as pl
import requests
import requests.adapters
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

try:
    from tqdm.notebook import tqdm as tqdm_notebook
except Exception:  # tqdm may be absent in some environments
    tqdm_notebook = None

from hoodini.models.schemas import GFF, NEIGHBORHOODS, PROTEINS, RECORDS
from hoodini.pipeline.helpers.neighborhood_extractor import extract_neighborhood
from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table
from hoodini.utils.logging_utils import error, info, success, warn
from hoodini.utils.polars_adapters import to_polars


def _extract_neighborhood_star(args):
    """Helper to allow imap_unordered with tuple args."""
    return extract_neighborhood(*args)


def in_jupyter():
    try:
        from IPython import get_ipython

        shell = get_ipython().__class__.__name__
        return shell == "ZMQInteractiveShell"
    except Exception:
        return False


def _enrich_proteins_with_metadata(
    all_prots: pl.DataFrame, records: pl.DataFrame, all_neigh: pl.DataFrame
) -> pl.DataFrame:
    """Enrich all_prots with target_prot (input protein ID), target_nuc (nucleotide ID), and uniprot_id.

    This ensures proteins have metadata linking them to their source protein and nucleotide,
    which is needed by downstream tools (padloc, defensefinder, cctyper, etc.) and visualization.

    Note: uniprot_id is only added to the actual target protein (where id matches protein_id),
    not to all proteins in the neighborhood.
    """
    if all_prots.is_empty():
        return all_prots

    if "unique_id" in all_prots.columns and "unique_id" in records.columns:
        # Join target_prot based on unique_id (applies to all proteins in neighborhood)
        prot_map = records.select(["unique_id", "protein_id"]).unique()
        prot_map = prot_map.with_columns(pl.col("unique_id").cast(pl.Utf8))
        all_prots = all_prots.with_columns(pl.col("unique_id").cast(pl.Utf8))

        all_prots = all_prots.join(
            prot_map.rename({"protein_id": "target_prot"}), on="unique_id", how="left"
        )

    # Join uniprot_id only to the actual target protein (where id == protein_id)
    if "uniprot_id" in records.columns and "id" in all_prots.columns:
        uniprot_map = (
            records.select(["protein_id", "uniprot_id"])
            .filter(pl.col("uniprot_id").is_not_null())
            .unique()
        )
        if uniprot_map.height > 0:
            all_prots = all_prots.join(uniprot_map, left_on="id", right_on="protein_id", how="left")

    if all_neigh is not None and all_neigh.height > 0:
        neigh_meta = all_neigh.select(["unique_id", "seqid"]).drop_nulls().unique()

        if "unique_id" in all_prots.columns and "unique_id" in neigh_meta.columns:
            neigh_meta = neigh_meta.with_columns(pl.col("unique_id").cast(pl.Utf8))

            all_prots = all_prots.join(
                neigh_meta.rename({"seqid": "target_nuc"}), on="unique_id", how="left"
            )

    return all_prots


def run_assembly_parser(
    records_df: pl.DataFrame,
    *,
    output_dir: Path | str | None = None,
    assembly_folder: str = None,
    ncrna: str | None = None,
    cctyper: bool = False,
    genomad: bool = False,
    blast: str = None,
    apikey: str = "",
    max_concurrent_downloads: int = 8,
    num_threads: int = 10,
    mod: str = "win_nts",
    wn: int = 20000,
    sorfs: bool = False,
    minwin: int = None,
    minwin_type: str = "both",
) -> dict:
    """
    Download assemblies and extract genomic neighborhoods around target proteins.

    Expected Files:
    ---------------
    - records_df: DataFrame from run_ipg with assembly accessions
    - assembly_folder: Optional pre-downloaded assemblies directory structure:
            {assembly_folder}/{GCA_XXXXXXXXX.X}/*.fna (genomic sequence)
            {assembly_folder}/{GCA_XXXXXXXXX.X}/*.gff (annotations)
            {assembly_folder}/{GCA_XXXXXXXXX.X}/*.faa (protein sequences)
    - Remote: NCBI FTP servers for assembly downloads if not in assembly_folder

    Generated Files:
    ----------------
    - {output}/assembly_list.txt: List of all assembly accessions to download
    - {output}/assembly_folder/{GCA_*}/*.fna: Downloaded genomic FASTA files
    - {output}/assembly_folder/{GCA_*}/*.gff: Downloaded GFF annotation files
    - {output}/assembly_folder/{GCA_*}/*.faa: Downloaded protein FASTA files
    - {output}/all_neigh.tsv: All extracted neighborhoods metadata
    - {output}/neighborhood/neighborhoods.fasta: Extracted neighborhood sequences
    - {output}/temp.gff: Temporary GFF of all extracted regions
    - {output}/results.fasta: All extracted protein sequences

    Process:
    --------
    1. Downloads assemblies (or uses local assembly_folder) for all valid assembly IDs
    2. Populates file paths (gbf_path, gff_path, faa_path, fna_path) for each record
    3. Runs extract_neighborhood() in parallel for each valid record
    4. Concatenates results into all_gff and all_neigh DataFrames
    5. Enforces minimum window size (minwin) and marks short contigs as failed
    6. Writes intermediate files for downstream stages

    Returns:
    --------
    dict with keys:
        - "records": updated DataFrame with file paths and extraction status
        - "all_gff": DataFrame of extracted GFF features (id, seqid, start, end, strand, etc.)
        - "all_prots": DataFrame of extracted protein sequences (id, sequence, product, etc.)
        - "all_neigh": DataFrame of neighborhood metadata (seqid, start_win, end_win, etc.)
        - "valid_uids": list of unique_id values for successfully extracted neighborhoods
    """
    records = to_polars(records_df, schema=RECORDS)

    if output_dir is None:
        input_path = None
        if "input_path" in records.columns:
            input_path = records.select("input_path").to_series()[0]
        if input_path is not None:
            output_dir = Path(input_path).with_suffix("")
        else:
            raise ValueError(
                "No output directory provided and could not infer from input_path. Please specify output_dir."
            )

    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve assembly folder: use provided path or default to output_dir/assembly_folder
    assembly_folder_path = (
        Path(assembly_folder).expanduser() if assembly_folder else output_dir / "assembly_folder"
    ).resolve()

    def _download_assemblies(df: pl.DataFrame) -> tuple[pl.DataFrame, str | None]:
        # Only download when we have an assembly_id, not failed, and we lack local files
        need_gbf = df["gbf_path"].is_null()
        need_struct = df["gff_path"].is_null() | df["faa_path"].is_null()  # fna optional
        mask = df["assembly_id"].is_not_null() & df["failed"].is_null() & need_gbf & need_struct
        assembly_list = df.filter(mask)["assembly_id"].drop_nulls().unique().to_list()
        # Normalize and keep only valid assembly accessions
        assembly_list = [
            str(aid).strip()
            for aid in assembly_list
            if isinstance(aid, str | Path) and str(aid).strip().startswith(("GCA_", "GCF_"))
        ]

        if not assembly_list:
            info("ℹ️  No assemblies to download.")
            return df, None

        # If using a provided assembly_folder, reuse already-downloaded GBFFs
        existing_map = {
            aid: str(assembly_folder_path / str(aid) / "genomic.gbff")
            for aid in assembly_list
            if (assembly_folder_path / str(aid) / "genomic.gbff").exists()
        }
        missing = [aid for aid in assembly_list if aid not in existing_map]
        if existing_map:
            df = df.with_columns(
                pl.when(pl.col("assembly_id").cast(pl.Utf8).is_in(list(existing_map.keys())))
                .then(pl.col("assembly_id").replace(existing_map))
                .otherwise(pl.col("gbf_path"))
                .alias("gbf_path")
            )
            info(
                f"✔️  Using {len(existing_map)} cached assemblies from {assembly_folder_path}"
                + (f" (missing: {missing})" if missing else "")
            )
            if not missing:
                return df, assembly_folder_path
            # keep only missing for download
            assembly_list = missing
        else:
            info("No cached assemblies found in assembly_folder; downloading all.")

        assembly_list_file = output_dir / "assembly_list.txt"
        with open(assembly_list_file, "w") as f:
            f.write("\n".join(assembly_list))

        links_result = get_prefetched_link_table(assembly_list)
        links = to_polars(links_result)

        if "filetype" in links.columns:
            gbff_links = links.filter(pl.col("filetype").str.to_lowercase() == "gbff")
        else:
            gbff_links = links

        missing_assemblies = set(assembly_list) - set(gbff_links["assembly_id"].unique())
        if missing_assemblies:
            warn(
                f"The following assembly_ids were not found in the download links: {missing_assemblies}"
            )

        session = requests.Session()
        total_files = gbff_links.height
        use_jupyter = in_jupyter()
        if use_jupyter:
            from tqdm.notebook import tqdm

            pbar = tqdm(total=total_files, desc="Downloading assemblies")
        else:
            from rich.progress import (
                BarColumn,
                Progress,
                SpinnerColumn,
                TaskProgressColumn,
                TextColumn,
                TimeElapsedColumn,
                TimeRemainingColumn,
            )

            progress = Progress(
                TextColumn(
                    f"[grey53][[/grey53][light_slate_grey]{datetime.now():%H:%M:%S}[/light_slate_grey][grey53]][/grey53]"
                ),
                SpinnerColumn(),
                TextColumn("{task.description}"),
                BarColumn(bar_width=40),
                TaskProgressColumn(),
                TextColumn("{task.completed}/{task.total}"),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
            )
            task = progress.add_task("Downloading assemblies", total=total_files)
            progress.start()

        def download_gbff(row: dict):
            assembly_id = row.get("assembly_id")
            url = row.get("url")
            target_folder = assembly_folder_path / str(assembly_id)
            target_file = target_folder / "genomic.gbff"

            if target_file.exists():
                return f"[SKIP] {assembly_id} already exists."

            try:
                target_folder.mkdir(parents=True, exist_ok=True)
                with session.get(url, timeout=60, stream=True) as resp:
                    resp.raise_for_status()
                    with open(target_file, "wb") as fout:
                        for chunk in resp.iter_content(chunk_size=1024 * 32):
                            if not chunk:
                                continue
                            fout.write(chunk)
                return f"[OK] {assembly_id}"
            except Exception as e:  # noqa: BLE001
                return f"[FAIL] {assembly_id}: {e}"

        max_threads = max(1, min(max_concurrent_downloads or 8, 32))
        adapter = requests.adapters.HTTPAdapter(
            pool_connections=max_threads, pool_maxsize=max_threads
        )
        session.mount("https://", adapter)
        session.mount("http://", adapter)

        link_records = gbff_links.to_dicts()
        with ThreadPoolExecutor(max_workers=max_threads) as executor:
            futures = [executor.submit(download_gbff, row) for row in link_records]
            for future in as_completed(futures):
                _ = future.result()
                if use_jupyter:
                    pbar.update(1)
                else:
                    progress.advance(task)
        if use_jupyter:
            pbar.close()
        else:
            progress.stop()

        # Polars replace expects primitive scalars; keep value as string path
        gbf_map = {
            aid: str(assembly_folder_path / str(aid) / "genomic.gbff") for aid in assembly_list
        }
        df = df.with_columns(
            pl.when(mask)
            .then(pl.col("assembly_id").replace(gbf_map))
            .otherwise(pl.col("gbf_path"))
            .alias("gbf_path")
        )

        success(f"Downloaded or located assemblies (folder: {assembly_folder_path})")
        return df, assembly_folder_path

    records, actual_assembly_folder = _download_assemblies(records)

    # Mark records without assembly_id and without local files as failed
    no_assembly_mask = (
        records["assembly_id"].is_null()
        & records["gbf_path"].is_null()
        & (records["gff_path"].is_null() | records["faa_path"].is_null())
        & records["failed"].is_null()
    )
    no_assembly_count = no_assembly_mask.sum()
    if no_assembly_count > 0:
        records = records.with_columns(
            pl.when(no_assembly_mask)
            .then(pl.lit(True))
            .otherwise(pl.col("failed"))
            .alias("failed"),
            pl.when(no_assembly_mask)
            .then(pl.lit("No assembly found for this record"))
            .otherwise(pl.col("failed_reason"))
            .alias("failed_reason"),
        )
        warn(f"{no_assembly_count} records have no assembly and were marked as failed")

    def write_fasta(df: pl.DataFrame, id_col: str, seq_col: str, path: str) -> None:
        with open(path, "w") as fh:
            for row in df.select([id_col, seq_col]).iter_rows(named=True):
                fh.write(f">{row[id_col]}\n")
                seq = row[seq_col]
                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")

    def _extract_neighborhoods(
        df: pl.DataFrame,
    ) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame, list]:
        df = df.with_row_count("unique_id").with_columns(pl.col("unique_id").cast(pl.Utf8))
        for col in ("start", "end", "strand"):
            if col not in df.columns:
                df = df.with_columns(pl.lit(None).alias(col))
        # Ensure failure flag/message columns exist for downstream joins
        if "failed" not in df.columns:
            df = df.with_columns(pl.lit(None).alias("failed"))
        if "failed_reason" not in df.columns:
            df = df.with_columns(pl.lit(None).alias("failed_reason"))

        file_list = []
        for row in df.iter_rows(named=True):
            has_gbf = bool(row.get("gbf_path"))
            has_gff_and_faa = bool(row.get("gff_path") and row.get("faa_path"))
            no_fail = row.get("failed") is None

            if (has_gbf or has_gff_and_faa) and no_fail:
                file_list.append(
                    (
                        row.get("protein_id"),
                        row.get("nucleotide_id"),
                        row.get("gbf_path"),
                        row.get("gff_path"),
                        row.get("faa_path"),
                        row.get("fna_path"),
                        mod,
                        wn,
                        row.get("strand"),
                        row.get("start"),
                        row.get("end"),
                        row.get("unique_id"),
                        row.get("input_type"),
                        sorfs,
                        row.get("is_full_contig", False),
                    )
                )

        if not file_list:
            info("ℹ️  No valid records to extract neighborhoods.")
            return df, pl.DataFrame(), pl.DataFrame(), pl.DataFrame(), []

        df_list: list[pl.DataFrame] = []
        df_list_neigh: list[pl.DataFrame] = []
        failed_ids: list[str] = []
        failed_msgs: list[str] = []

        use_jupyter = in_jupyter()
        pbar = None
        progress = None
        task = None
        try:
            if use_jupyter:
                if tqdm_notebook is not None:
                    pbar = tqdm_notebook(total=len(file_list), desc="Parsing GBFF")
                else:
                    warn("tqdm not available; falling back to Rich progress.")
                    use_jupyter = False
            else:
                progress = Progress(
                    TextColumn(
                        f"[grey53][[/grey53][light_slate_grey]{datetime.now():%H:%M:%S}[/light_slate_grey][grey53]][/grey53]"
                    ),
                    SpinnerColumn(),
                    TextColumn("{task.description}"),
                    BarColumn(bar_width=40),
                    TaskProgressColumn(),
                    TextColumn("{task.completed}/{task.total}"),
                    TimeElapsedColumn(),
                    TimeRemainingColumn(),
                )
                task = progress.add_task("Parsing GBFF", total=len(file_list))
                progress.start()

            with ProcessPoolExecutor(
                max_workers=num_threads, mp_context=_mp.get_context("spawn")
            ) as executor:
                futures = {
                    executor.submit(_extract_neighborhood_star, item): idx
                    for idx, item in enumerate(file_list)
                }

                for completed, future in enumerate(as_completed(futures), start=1):
                    idx = futures[future]
                    item = file_list[idx]

                    try:
                        res = future.result(timeout=600)
                        if res[0] is not None:
                            df_list.append(to_polars(res[0]))
                        if res[1] is not None:
                            df_list_neigh.append(to_polars(res[1]))
                        if res[0] is None:
                            failed_ids.append(res[2])
                            failed_msgs.append(res[3])
                    except TimeoutError:
                        failed_ids.append(str(item))
                        failed_msgs.append("TIMEOUT after 10 minutes")
                    except Exception as e:
                        failed_ids.append(str(item))
                        failed_msgs.append(f"Exception: {type(e).__name__}: {str(e)}")
                    finally:
                        if use_jupyter and pbar:
                            pbar.update(1)
                        elif progress and task is not None:
                            progress.advance(task)
        finally:
            if pbar:
                pbar.close()
            if progress:
                progress.stop()

        all_prots = pl.concat(df_list, how="vertical") if df_list else pl.DataFrame()
        if not all_prots.is_empty() and "id" in all_prots.columns:
            all_prots = all_prots.with_columns((pl.lit("ID=") + pl.col("id")).alias("attributes"))

        gff_columns = [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
        all_gff = (
            all_prots.select([c for c in gff_columns if c in all_prots.columns])
            if not all_prots.is_empty()
            else pl.DataFrame()
        )
        all_gff = to_polars(all_gff, schema=GFF) if not all_gff.is_empty() else all_gff
        if not all_prots.is_empty():
            all_prots = all_prots.drop([c for c in gff_columns if c in all_prots.columns])
            all_prots = to_polars(all_prots, schema=PROTEINS)

        all_neigh = pl.concat(df_list_neigh, how="vertical") if df_list_neigh else pl.DataFrame()
        all_neigh = (
            to_polars(all_neigh, schema=NEIGHBORHOODS) if not all_neigh.is_empty() else all_neigh
        )

        if not all_neigh.is_empty():
            all_neigh = all_neigh.with_columns(
                (pl.col("end_win") - pl.col("start_win")).alias("length")
            )
            all_neigh.drop("sequence").write_csv(output_dir / "all_neigh.tsv", separator="\t")

            if minwin is not None:
                # Exclude full contig analyses from minwin filtering
                # (when start/end weren't specified, the whole contig is the target)
                if "is_full_contig" in all_neigh.columns:
                    filterable = all_neigh.filter(pl.col("is_full_contig") == False)  # noqa: E712
                else:
                    filterable = all_neigh
                
                if minwin_type == "total":
                    short_contigs = (
                        filterable.filter(pl.col("length") < minwin)["unique_id"].unique().to_list()
                    )
                elif minwin_type == "upstream":
                    short_contigs = (
                        filterable.filter((pl.col("start_target") - pl.col("start_win")) < minwin)[
                            "unique_id"
                        ]
                        .unique()
                        .to_list()
                    )
                elif minwin_type == "downstream":
                    short_contigs = (
                        filterable.filter((pl.col("end_win") - pl.col("end_target")) < minwin)[
                            "unique_id"
                        ]
                        .unique()
                        .to_list()
                    )
                else:
                    short_contigs = (
                        filterable.filter(
                            ((pl.col("start_target") - pl.col("start_win")) < minwin)
                            | ((pl.col("end_win") - pl.col("end_target")) < minwin)
                        )["unique_id"]
                        .unique()
                        .to_list()
                    )
            else:
                short_contigs = []

            neigh_dir = output_dir / "neighborhood"
            neigh_dir.mkdir(parents=True, exist_ok=True)

            if failed_ids:
                warn(f"Failed extractions for {len(failed_ids)} records: {failed_ids}")

            if short_contigs:
                df = (
                    df.join(
                        pl.DataFrame(
                            {
                                "unique_id": short_contigs,
                                "failed": [True] * len(short_contigs),
                                "failed_reason": [
                                    "Genomic context shorter than minimum window size"
                                ]
                                * len(short_contigs),
                            }
                        ),
                        on="unique_id",
                        how="left",
                        suffix="_short",
                    )
                    .with_columns(
                        pl.when(pl.col("failed_short").is_not_null())
                        .then(pl.col("failed_short"))
                        .otherwise(pl.col("failed"))
                        .alias("failed"),
                        pl.when(pl.col("failed_reason_short").is_not_null())
                        .then(pl.col("failed_reason_short"))
                        .otherwise(pl.col("failed_reason"))
                        .alias("failed_reason"),
                    )
                    .drop("failed_short", "failed_reason_short")
                )
                warn(f"{len(short_contigs)} contigs below min window size: {short_contigs}")

            # Check if any valid records remain after marking failures
            valid_count = df.filter(pl.col("failed").is_null()).height
            total_count = df.height
            failed_count = total_count - valid_count
            # Count pre-existing failures (e.g., no assembly)
            pre_existing_failures = failed_count - len(short_contigs) - len(failed_ids)
            if valid_count == 0:
                df.write_csv(output_dir / "records.tsv", separator="\t", quote_style="necessary")
                failure_parts = []
                if len(short_contigs) > 0:
                    failure_parts.append(f"{len(short_contigs)} below min window size")
                if len(failed_ids) > 0:
                    failure_parts.append(f"{len(failed_ids)} extraction failures")
                if pre_existing_failures > 0:
                    failure_parts.append(f"{pre_existing_failures} without assembly")
                failure_summary = ", ".join(failure_parts) if failure_parts else "unknown reasons"
                error(
                    f"Aborting! All {total_count} records failed ({failure_summary}). "
                    f"Consider decreasing --min-win to accept smaller neighborhoods, or checking your input sequences."
                )
                sys.exit(1)

            if failed_ids:
                df = (
                    df.join(
                        pl.DataFrame(
                            {
                                "unique_id": failed_ids,
                                "failed": [True] * len(failed_ids),
                                "failed_reason": failed_msgs,
                            }
                        ),
                        on="unique_id",
                        how="left",
                        suffix="_fail",
                    )
                    .with_columns(
                        pl.when(pl.col("failed_fail").is_not_null())
                        .then(pl.col("failed_fail"))
                        .otherwise(pl.col("failed"))
                        .alias("failed"),
                        pl.when(pl.col("failed_reason_fail").is_not_null())
                        .then(pl.col("failed_reason_fail"))
                        .otherwise(pl.col("failed_reason"))
                        .alias("failed_reason"),
                    )
                    .drop("failed_fail", "failed_reason_fail")
                )

            valid_uids = (
                df.filter(pl.col("nucleotide_id").is_not_null() & pl.col("failed").is_null())[
                    "unique_id"
                ]
                .unique()
                .to_list()
            )

            all_neigh = all_neigh.with_columns(
                (
                    pl.col("seqid").cast(pl.Utf8)
                    + pl.lit("_")
                    + pl.col("start_win").cast(pl.Utf8)
                    + pl.lit("_")
                    + pl.col("end_win").cast(pl.Utf8)
                ).alias("temp_seqid")
            )
            subset_neigh = (
                all_neigh.filter(pl.col("unique_id").is_in(valid_uids))
                .select(["temp_seqid", "sequence"])
                .drop_nulls()
                .unique(subset=["temp_seqid"])
            )

            neigh_fasta = neigh_dir / "neighborhoods.fasta"
            if subset_neigh.height > 0:
                write_fasta(subset_neigh, "temp_seqid", "sequence", neigh_fasta)

            df.write_csv(output_dir / "records.tsv", separator="\t", quote_style="necessary")

            subset_gff = (
                all_prots.select([c for c in ["id", "sequence"] if c in all_prots.columns])
                .drop_nulls()
                .unique(subset=["id"])
            )
            results_fasta = output_dir / "results.fasta"
            if subset_gff.height > 0:
                write_fasta(subset_gff, "id", "sequence", results_fasta)

            success("Extracted neighborhoods and wrote output files.")
            return df, all_gff, all_prots, all_neigh, valid_uids

        else:
            if failed_ids:
                df = (
                    df.join(
                        pl.DataFrame(
                            {
                                "unique_id": failed_ids,
                                "failed": [True] * len(failed_ids),
                                "failed_reason": failed_msgs,
                            }
                        ),
                        on="unique_id",
                        how="left",
                        suffix="_fail",
                    )
                    .with_columns(
                        pl.when(pl.col("failed_fail").is_not_null())
                        .then(pl.col("failed_fail"))
                        .otherwise(pl.col("failed"))
                        .alias("failed"),
                        pl.when(pl.col("failed_reason_fail").is_not_null())
                        .then(pl.col("failed_reason_fail"))
                        .otherwise(pl.col("failed_reason"))
                        .alias("failed_reason"),
                    )
                    .drop("failed_fail", "failed_reason_fail")
                )
            df.write_csv(output_dir / "records.tsv", separator="\t", quote_style="necessary")
            error("Aborting! All IDs failed to yield neighborhoods.")
            sys.exit(1)

    updated_records, all_gff, all_prots, all_neigh, valid_uids = _extract_neighborhoods(records)

    all_prots = _enrich_proteins_with_metadata(all_prots, updated_records, all_neigh)

    if (
        not all_prots.is_empty()
        and "protein_id" in all_prots.columns
        and "target_prot" in all_prots.columns
    ):
        add_data = all_prots.filter(pl.col("protein_id") == pl.col("target_prot")).select(
            ["unique_id", "target_prot"]
        )
        all_neigh = all_neigh.join(
            add_data.with_columns(pl.col("unique_id").cast(pl.Utf8)), on="unique_id", how="left"
        )

    return {
        "records": updated_records,
        "all_gff": all_gff,
        "all_prots": all_prots,
        "all_neigh": all_neigh,
        "valid_uids": valid_uids,
    }
