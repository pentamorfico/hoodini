import gc
import subprocess
from importlib.resources import files
from pathlib import Path
from shutil import copyfile

import polars as pl
import pyarrow.parquet as pq

from hoodini.utils.logging_utils import info, success, warn


def _stream_filter_parquet_by_id(
    path: str, ids: set[int], columns: list[str], id_col: str = "id"
) -> pl.DataFrame:
    """Filter a large parquet file by id, one row group at a time to bound peak memory."""
    pf = pq.ParquetFile(path)
    chunks: list[pl.DataFrame] = []
    for rg_idx in range(pf.num_row_groups):
        table = pf.read_row_group(rg_idx, columns=columns)
        df = pl.from_arrow(table)
        del table
        filtered = df.filter(pl.col(id_col).is_in(ids))
        del df
        if filtered.height:
            # compact copy; matched rows must not pin the parent row-group buffer
            chunks.append(pl.from_arrow(filtered.to_arrow().combine_chunks()))
        del filtered
        gc.collect()

    if not chunks:
        return pl.DataFrame(schema=dict.fromkeys(columns, pl.Null))
    return pl.concat(chunks)


def run_emapper(all_prots: pl.DataFrame, output: str | Path, num_threads: int = 1) -> pl.DataFrame:
    """
    Run DIAMOND blastp, pick best hit per query directly in Polars,
    join to eggNOG metadata, pick the deepest OG per query,
    and return one row per input protein as a Polars DataFrame.

    Streams the eggnog_prots.parquet lookup via `_stream_filter_parquet_by_id`.
    """

    info("🧾\tRunning eggNOG-mapper (DIAMOND + eggNOG, best+deepest OG) ...")

    output = Path(output)
    emapper_dir = output / "emapper"
    emapper_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = output / "results.faa"
    fasta_fallback = output / "results.fasta"

    if not fasta_path.exists():
        if fasta_fallback.exists():
            copyfile(fasta_fallback, fasta_path)
            info(f"Copied {fasta_fallback} -> {fasta_path}")
        else:
            seq_df = all_prots[["id", "sequence"]].drop_nulls().drop_duplicates("id")
            seq_df.to_fasta("id", "sequence", fasta_path)
            success(f"Generated {fasta_path}")

    diamond_db = str(files("hoodini").joinpath("data", "emapper", "eggnog_proteins.dmnd"))

    results_m8 = emapper_dir / "results.m8"

    cmd = [
        "diamond",
        "blastp",
        "-q",
        str(fasta_path),
        "-d",
        diamond_db,
        "-o",
        str(results_m8),
        "--threads",
        str(max(1, int(num_threads or 1))),
        "--max-target-seqs",
        "1",
        "--evalue",
        "0.001",
    ]
    info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    if not results_m8.exists():
        warn(f"DIAMOND results not found at {results_m8}")
        return pl.DataFrame()

    hits_all = pl.read_csv(
        results_m8,
        has_header=False,
        separator="\t",
        new_columns=[
            "qseqid",
            "sseqid",
            "pident",
            "alnlen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )

    hits_best = (
        hits_all.sort(["qseqid", "bitscore"], descending=[False, True])
        .group_by("qseqid")
        .head(1)
        .select(["qseqid", "sseqid"])
        .with_columns(pl.col("sseqid").cast(pl.Int64))
    )

    eggnog_prots_path = str(files("hoodini").joinpath("data", "emapper", "eggnog_prots.parquet"))
    eggnog_og_path = str(files("hoodini").joinpath("data", "emapper", "eggnog_og.parquet"))

    # Get list of sseqids we need to look up. DIAMOND/MMseqs2 report the
    # subject as eggNOG's internal numeric protein id (eggnog_prots.parquet's
    # `id` column), not the human-readable "taxid.locus" `name` string, so the
    # lookup below joins on `id`.
    prot_ids = hits_best["sseqid"].unique().to_list()

    # eggNOG 7: eggnog_prots.parquet's `ogs` column is a comma-separated list of
    # full OG identifiers ("cluster@taxid|clade[!]"), each of which matches
    # eggnog_og.parquet's `name` column exactly, so we join on `name` directly
    # instead of splitting into separate (og, level) parts and joining on both.
    prot_cols = [
        "gos",
        "pfam",
        "kegg_ko",
        "kegg_ec",
        "kegg_pathway",
        "kegg_module",
        "kegg_reaction",
        "kegg_rclass",
        "kegg_brite",
        "kegg_tc",
        "kegg_cazy",
        "kegg_cog",
        "kegg_disease",
        "kegg_go",
        "kegg_drug",
        "kegg_pubmed",
        "kegg_network",
        "bigg_reaction",
    ]

    # Streaming lookup; see _stream_filter_parquet_by_id for why not DuckDB.
    prot_id_set = set(prot_ids)
    prots = _stream_filter_parquet_by_id(
        eggnog_prots_path,
        prot_id_set,
        columns=["id", "name", "ogs", "pname"] + prot_cols,
    ).rename({"id": "prot_id"})

    prots = (
        prots.with_columns(pl.col("ogs").fill_null("").str.split(",").alias("og_name"))
        .explode("og_name")
        .filter(pl.col("og_name") != "")
        .drop("ogs")
    )
    og = pl.read_parquet(eggnog_og_path).rename({"pname": "og_pname"})
    annotated = prots.join(og, left_on="og_name", right_on="name", how="inner").drop("og_name")

    hits_annotated = hits_best.join(annotated, left_on="sseqid", right_on="prot_id", how="left")

    # Pick the deepest OG (highest `depth`) per query protein.
    deepest = (
        hits_annotated.sort(["qseqid", "depth"], descending=[False, True], nulls_last=True)
        .group_by("qseqid")
        .head(1)
    )

    deepest = deepest.rename({"qseqid": "id"})
    exclude_cols = {"level", "nm", "sseqid", "name"}
    lead = ["id", "pname", "description", "COG_categories", "pfam"]
    lead_present = [c for c in lead if c in deepest.columns]
    rest = [c for c in deepest.columns if c not in lead_present and c not in exclude_cols]
    deepest = deepest.select(lead_present + rest)

    info("🔎 Head of annotated Polars DF (deepest OG per best hit):")
    info(deepest.head(10))
    info(f"shape: {deepest.shape}")

    success(f"DIAMOND annotations ready: {deepest.height} queries annotated")
    return deepest
