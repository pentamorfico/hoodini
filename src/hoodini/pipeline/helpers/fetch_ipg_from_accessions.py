import argparse
import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from io import StringIO
from pathlib import Path

import polars as pl
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from hoodini.utils.logging_utils import error, info, warn
from hoodini.utils.ncbi_eutils import efetch_with_retries

TOOL = "ipg_fetcher"
CHUNK_SIZE = 100

DEFAULT_MAX_CONCURRENT = 9
DEFAULT_FALLBACK_CONCURRENT = 3
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
MAX_WORKERS = DEFAULT_MAX_CONCURRENT if NCBI_API_KEY else DEFAULT_FALLBACK_CONCURRENT
MAX_PARALLEL = MAX_WORKERS


def _efetch_chunk(accessions: list[str]) -> str:
    try:
        return efetch_with_retries(
            db="protein",
            ids=accessions,
            rettype="ipg",
            retmode="text",
            api_key=NCBI_API_KEY,
            tool=TOOL,
            timeout=90,
            max_retries=3,
            backoff_seconds=5,
        )
    except RuntimeError as e:
        error(f"efetch failed for IDs {accessions[:3]}...: {e}")
        return ""


def fetch_ipg_from_accessions(accessions: list[str]) -> pl.DataFrame:
    """
    Fetch IPG data for a list of protein accessions in parallel with rate limiting.
    """
    semaphore = threading.Semaphore(MAX_PARALLEL)
    chunks = [accessions[i : i + CHUNK_SIZE] for i in range(0, len(accessions), CHUNK_SIZE)]
    results = [None] * len(chunks)

    def wrapped_efetch(chunk, idx):
        with semaphore:
            return idx, _efetch_chunk(chunk)

    ts_col = TextColumn(
        f"[grey53][[/grey53][light_slate_grey]{datetime.now():%H:%M:%S}[/light_slate_grey][grey53]][/grey53]"
    )
    with Progress(
        ts_col,
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(bar_width=40),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("Chunk {task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task("Fetching IPG chunks", total=len(chunks))
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            future_to_idx = {
                executor.submit(wrapped_efetch, chunk, idx): idx for idx, chunk in enumerate(chunks)
            }
            for future in as_completed(future_to_idx):
                idx, chunk_result = future.result()
                if chunk_result and chunk_result.strip():
                    try:
                        chunk_df = pl.read_csv(
                            StringIO(chunk_result),
                            separator="\t",
                            infer_schema_length=1000,
                        )
                        chunk_df = chunk_df.rename({col: col.lower() for col in chunk_df.columns})
                        results[idx] = chunk_df
                    except Exception as e:
                        snippet = chunk_result[:200].replace("\n", "\\n")
                        warn(f"Failed to parse IPG chunk {idx}: {e}. Snippet: {snippet}")
                        results[idx] = None
                else:
                    warn(f"Empty IPG chunk {idx} returned by efetch")
                    results[idx] = None
                progress.update(task, advance=1)

    dfs = [df for df in results if df is not None and df.height > 0]
    if not dfs:
        return pl.DataFrame()

    return pl.concat(dfs, how="vertical_relaxed")


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve IPG data via NCBI's efetch.fcgi for protein accessions."
    )
    parser.add_argument("input", help="Input file with protein accessions (one per line)")
    parser.add_argument("-o", "--output", help="Output TSV file to save IPG results", required=True)
    args = parser.parse_args()

    with open(args.input) as f:
        accessions = [line.strip() for line in f if line.strip()]

    df = fetch_ipg_from_accessions(accessions)

    if df.height == 0:
        warn("No IPG data retrieved.")
    else:
        out_path = Path(args.output)
        df.write_csv(out_path, separator="\t")
        info(f"Saved {df.height} IPG entries to {out_path}")


if __name__ == "__main__":
    main()
