import datetime
import shutil
import sys
from importlib.resources import files
from pathlib import Path

import polars as pl

from hoodini.models.schemas import RECORDS
from hoodini.pipeline.helpers.single_query import prepare_single_query_input
from hoodini.utils.logging_utils import error, info, warn
from hoodini.utils.polars_adapters import to_polars
from hoodini.utils.validation import read_input_list, read_input_sheet, uniprot2ncbi


def initialize_inputs(
    *,
    input_path: Path | str | None = None,
    inputsheet: Path | str | None = None,
    output: Path | str | None = None,
    force: bool = False,
    remote_evalue: float = 1e-5,
    remote_max_targets: int = 100,
) -> pl.DataFrame:
    """
    Initialize the working directory and read the user's input records (Polars).

    Expected Files:
    ---------------
    - input_path: single-column text file with protein IDs (NCBI/UniProt), one per line
        OR
    - inputsheet: TSV file with columns: og_index, seqid, accession, organism, etc.
    - hoodini/data/assembly_summary.parquet (packaged with hoodini, auto-checked)

    Generated Files:
    ----------------
    - {output}/ directory (created or overwritten if force=True)
    - No immediate output files; returns DataFrame for downstream stages

    Process:
    --------
    1. Creates or (if existing) optionally overwrites the output folder.
    2. Reads either a single‐column input list or a TSV "inputsheet".
    3. Converts UniProt IDs to NCBI IDs via `uniprot2ncbi(...)` with remote BLAST if needed.
    4. Drops duplicate records based on "og_index".
    5. Returns a Polars DataFrame of final, deduplicated records.

    Returns:
    --------
    pl.DataFrame with schema matching RECORDS (og_index, seqid, accession, organism, etc.)
    """

    check_assembly_db()
    check_contig_lengths_db()

    if output:
        output_folder = Path(output)
    elif input_path:
        output_folder = Path(input_path).with_suffix("")
    elif inputsheet:
        output_folder = Path(inputsheet).with_suffix("")
    else:
        error("Either `input_path` or `inputsheet` must be provided.")
        sys.exit(1)

    if output_folder.exists():
        if force:
            warn(f"Overwriting existing folder {output_folder}.")
            shutil.rmtree(output_folder)
        else:
            error(f"Output folder '{output_folder}' already exists. Use --force to overwrite.")
            sys.exit(1)

    output_folder.mkdir(parents=True, exist_ok=True)
    info(f"✔️  Created folder {output_folder}")

    # If input_path is a literal (not a file), build a temp list inside the output folder.
    if input_path and not Path(input_path).exists():
        temp_input = prepare_single_query_input(
            str(input_path),
            output_folder,
            evalue=remote_evalue,
            max_targets=remote_max_targets,
        )
        if not temp_input:
            error("Failed to prepare input from single query; aborting.")
            sys.exit(1)
        input_path = str(temp_input)

    if inputsheet:
        records_raw = read_input_sheet(inputsheet)
    elif input_path:
        records_raw = read_input_list(input_path)
    else:
        error("No input_path or inputsheet provided.")
        sys.exit(1)

    records_pd = uniprot2ncbi(records_raw)
    records = to_polars(records_pd, schema=RECORDS)

    records = records.unique(subset=["og_index"], keep="first")

    records = records.with_columns(
        (
            (pl.col("gff_path").is_not_null() & pl.col("faa_path").is_not_null())
            | (pl.col("gbf_path").is_not_null())
        ).alias("premade")
    )

    return records


def check_assembly_db() -> None:
    """
    Check if assembly_summary.parquet exists and optionally download it.
    If the file is older than 1 month, show a warning and instruct the user to run the update command.
    """
    try:
        summary_path = files("hoodini").joinpath("data", "assembly_summary.parquet")
        now = datetime.datetime.now()
        one_month = datetime.timedelta(days=30)

        if summary_path.exists():
            mtime = datetime.datetime.fromtimestamp(summary_path.stat().st_mtime)
            age = now - mtime
            info(f"📁 Assembly DB found: {summary_path} (updated: {mtime:%Y-%m-%d})")
            if age > one_month:
                warn(
                    "The assembly database is older than 30 days. Run 'hoodini download assembly_summary' to update."
                )
            else:
                info("✅ Using existing database (less than 1 month old).")
            return
        else:
            info("📭 No local NCBI assembly summary database found. Downloading now...")
            from hoodini.download.assembly_summary import download_assembly_summary_db

            output_path = download_assembly_summary_db(summary_path)
            info(f"✔️  Downloaded assembly database to {output_path}")
            return
    except Exception as e:
        error(f"Error checking or downloading assembly DB: {e}")


def check_contig_lengths_db() -> None:
    """
    Check if contig_lengths parquet files exist and download if missing.
    Downloads from remote storage if no local files are found.
    """
    REMOTE_URL = "https://storage.hoodini.bio/contig_lengths.parquet"

    try:
        contig_dir = files("hoodini").joinpath("data", "contig_lengths")

        # Check if any parquet files exist in the directory
        parquet_files = list(contig_dir.glob("*.parquet")) if contig_dir.exists() else []

        if parquet_files:
            info(
                f"📁 Contig lengths DB found: {len(parquet_files)} parquet file(s) in {contig_dir}"
            )
            return

        # No files found - download from remote
        info("📭 No local contig lengths database found. Downloading now...")
        info("   This may take a few minutes (downloading ~2GB file)...")

        from hoodini.download.databases import _download_url

        # Ensure directory exists
        contig_dir.mkdir(parents=True, exist_ok=True)

        # Download to contig_lengths.parquet
        dest = contig_dir / "contig_lengths.parquet"
        success = _download_url(REMOTE_URL, dest)

        if success and dest.exists():
            info(f"✔️  Downloaded contig lengths database to {dest}")
        else:
            warn(
                "⚠️  Failed to download contig_lengths.parquet. "
                "Run 'hoodini download contig_lengths' manually to retry."
            )

    except Exception as e:
        error(f"Error checking or downloading contig lengths DB: {e}")
