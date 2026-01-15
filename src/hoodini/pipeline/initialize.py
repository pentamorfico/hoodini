import datetime
import shutil
import sys
from importlib.resources import files
from pathlib import Path

import polars as pl

from hoodini.models.schemas import RECORDS
from hoodini.pipeline.helpers.single_query import prepare_single_query_input
from hoodini.utils.logging_utils import error, info, prompt, warn
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
    check_playwright_browser()

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
            warn(f"Folder {output_folder} already exists.")
            answer = _prompt_yes_no("⌨️  Remove it?", default="y/N")
            if not answer:
                error("Aborting: Folder not modified.")
                sys.exit(1)
            shutil.rmtree(output_folder)

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


def check_playwright_browser() -> None:
    """
    Check if Playwright Chromium is installed and install it if missing.
    """
    import os
    import subprocess
    import sys
    from pathlib import Path

    # Skip in Docker or CI
    if Path("/.dockerenv").exists() or os.getenv("CI"):
        return

    try:
        # Standard playwright browser cache locations
        home = Path.home()
        possible_paths = [
            home / ".cache" / "ms-playwright",  # Linux
            home / "Library" / "Caches" / "ms-playwright",  # macOS
            home / "AppData" / "Local" / "ms-playwright",  # Windows
        ]
        
        # Check if chromium exists in any cache location
        for cache_path in possible_paths:
            if cache_path.exists():
                chromium_dirs = list(cache_path.glob("chromium*"))
                if chromium_dirs:
                    return  # Already installed
        
        info("🎭 Playwright Chromium not found. Installing (one-time setup)...")
        subprocess.run(
            [sys.executable, "-m", "playwright", "install", "--only-shell", "chromium"],
            check=True,
            timeout=180,
        )
        info("✅ Playwright Chromium installed successfully")
            
    except subprocess.TimeoutExpired:
        warn("Playwright install timed out. Run 'playwright install chromium' manually.")
    except Exception as e:
        warn(f"Could not check/install Playwright: {e}")



def _prompt_yes_no(prompt_text: str, default: str = "no") -> bool:
    """Ask user a yes/no question via prompt."""
    valid_yes = {"y", "yes"}
    valid_no = {"n", "no"}
    while True:
        response = prompt(prompt_text, default=default).strip().lower()
        if response in valid_yes:
            return True
        if response in valid_no:
            return False
        warn("Please enter 'y' or 'n'.")
