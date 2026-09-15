import asyncio
import contextlib
import os
from datetime import datetime
from email.utils import parsedate_to_datetime
from importlib.resources import files
from pathlib import Path

import duckdb
import polars as pl  # type: ignore[import]
import pyarrow.parquet as pq  # type: ignore[import]
import requests  # type: ignore[import]

from hoodini.download.assembly_summary import download_assembly_summary_db
from hoodini.download.ncbi_sequence_reports import PartRotatingWriter, fetch_sequence_reports
from hoodini.utils.contig_metadata import register_contig_table, scan_contig_table
from hoodini.utils.logging_utils import console

NCBI_API_KEY: str | None = os.environ.get("NCBI_API_KEY")

DATA_DIR = files("hoodini").joinpath("data")
CONTIG_LENGTHS_DIR = DATA_DIR.joinpath("contig_lengths")
MASTER_CONTIGS = DATA_DIR.joinpath("contig_lengths.parquet")
ASSEMBLY_SUMMARY = DATA_DIR.joinpath("assembly_summary.parquet")
CONTIG_LENGTHS_CHECKPOINT = DATA_DIR.joinpath("contig_lengths_checkpoint.sqlite")

DEFAULT_GROUPS = {"bacteria", "viral", "archaea", "metagenomes", "other"}
MAX_RETRIES = 5
REMOTE_CONTIG_LENGTHS_URL = "https://storage.hoodini.bio/contig_lengths.parquet"

__all__ = [
    "PartRotatingWriter",
    "get_missing_contigs_from_summary",
    "download_contig_lengths",
]

_ASM_CANDIDATES: tuple[str, ...] = (
    "assembly_accession",
    "assemblyAccession",
    "assembly_id",
    "assemblyId",
    "assembly",
)


def _get_remote_parquet_last_modified() -> datetime | None:
    """Get Last-Modified date from remote contig_lengths.parquet."""
    try:
        resp = requests.head(REMOTE_CONTIG_LENGTHS_URL, timeout=10)
        if resp.status_code == 200:
            last_mod = resp.headers.get("Last-Modified")
            if last_mod:
                dt = parsedate_to_datetime(last_mod)
                console.log(f"Remote contig_lengths.parquet last modified: {dt}")
                return dt
    except Exception as e:
        console.log(f"⚠️  Could not fetch remote Last-Modified: {e}")
    return None


def _detect_assembly_col_in_file(path: Path) -> str | None:
    try:
        pf = pq.ParquetFile(path)
        names = set(pf.schema.names)
        for c in _ASM_CANDIDATES:
            if c in names:
                return c
    except Exception:
        pass
    return None


def _detect_assembly_col_in_dir(dirpath: Path, sample: int = 10) -> str | None:
    try:
        files = sorted(dirpath.glob("part-*.parquet"))[:sample]
        for f in files:
            c = _detect_assembly_col_in_file(f)
            if c is not None:
                return c
    except Exception:
        pass
    return None


def get_missing_contigs_from_summary(
    assembly_summary_path: Path,
    allowed_assemblies_df: pl.DataFrame | None = None,
) -> tuple[pl.DataFrame, float | None]:
    """Return a DataFrame of missing assembly_accession values and latest mtime.

    Uses DuckDB for memory-efficient querying of large parquet files.
    """
    latest_mtime: float | None = None
    all_parquet_files = list(CONTIG_LENGTHS_DIR.glob("*.parquet"))

    if all_parquet_files:
        console.log(f"Scanning {len(all_parquet_files)} parquet files in contig_lengths/...")
        with contextlib.suppress(Exception):
            latest_mtime = max(p.stat().st_mtime for p in all_parquet_files)

    try:
        with duckdb.connect(":memory:") as con:
            con.execute('SET memory_limit = "4GB"')

            # Expose allowed assemblies as a zero-copy DuckDB view. A prior
            # version built this via CREATE TABLE + row-by-row executemany,
            # which takes minutes for the several-million-row candidate lists
            # this receives; register() scans the Polars DataFrame directly.
            if allowed_assemblies_df is not None:
                con.register("allowed_asm", allowed_assemblies_df)

            # Build the query for valid assemblies from assembly_summary
            groups_str = ", ".join(f"'{g}'" for g in DEFAULT_GROUPS)
            summary_query = f"""
                SELECT DISTINCT CAST(assembly_accession AS VARCHAR) as assembly_accession
                FROM read_parquet('{str(assembly_summary_path)}')
                WHERE "group" IN ({groups_str})
                  AND ftp_path IS NOT NULL
                  AND TRIM(ftp_path) != ''
                  AND LOWER(ftp_path) != 'na'
            """

            if allowed_assemblies_df is not None:
                summary_query += (
                    " AND assembly_accession IN (SELECT assembly_accession FROM allowed_asm)"
                )

            if all_parquet_files:
                # Query existing contig assemblies and do anti-join
                contig_glob = str(CONTIG_LENGTHS_DIR / "*.parquet")
                register_contig_table(con, contig_glob)
                missing_df = con.execute(
                    f"""
                    WITH summary AS ({summary_query}),
                    existing AS (
                        SELECT DISTINCT CAST(assemblyAccession AS VARCHAR) as assembly_accession
                        FROM hoodini_contigs
                    )
                    SELECT s.assembly_accession
                    FROM summary s
                    LEFT JOIN existing e ON s.assembly_accession = e.assembly_accession
                    WHERE e.assembly_accession IS NULL
                """
                ).pl()
            else:
                console.log("No existing contig_lengths found, will download all")
                missing_df = con.execute(summary_query).pl()

    except Exception as e:
        console.log(f"⚠️  DuckDB failed, falling back to Polars streaming: {e}")
        # Fallback to Polars if DuckDB fails
        summary_lf = (
            pl.scan_parquet(str(assembly_summary_path))
            .filter(
                (pl.col("group").is_in(DEFAULT_GROUPS))
                & pl.col("ftp_path").is_not_null()
                & (pl.col("ftp_path").str.strip_chars() != "")
                & (pl.col("ftp_path").str.to_lowercase() != "na")
            )
            .select(pl.col("assembly_accession").cast(pl.Utf8))
            .unique()
        )
        if allowed_assemblies_df is not None:
            summary_lf = summary_lf.join(
                allowed_assemblies_df.lazy(),
                on="assembly_accession",
                how="semi",
            )
        if all_parquet_files:
            contig_lf = (
                scan_contig_table(
                    str(CONTIG_LENGTHS_DIR / "*.parquet"),
                )
                .select(pl.col("assemblyAccession").cast(pl.Utf8).alias("assembly_accession"))
                .unique()
            )
            missing_lf = summary_lf.join(contig_lf, on="assembly_accession", how="anti")
        else:
            missing_lf = summary_lf
        missing_df = missing_lf.collect(streaming=True)

    console.log(f"❗ {missing_df.height:,} assemblies missing contig lengths.")
    console.log(f"latest_mtime: {latest_mtime}")

    return missing_df, latest_mtime


def _validate_api_keys(keys: list[str]) -> list[str]:
    """Probe each NCBI API key with a cheap request and drop any that fail.

    Prevents malformed/typo'd/revoked keys (which fail every request) from
    silently eating retries and reducing effective throughput once assigned
    to a rotation slot.
    """
    from hoodini.download.ncbi_sequence_reports import SEQUENCE_REPORTS_URL

    valid = []
    for key in keys:
        suffix = key[-4:] if len(key) >= 4 else key
        try:
            resp = requests.post(
                SEQUENCE_REPORTS_URL,
                json={"accession": "GCF_000005845.2", "page_size": 1},
                headers={"api-key": key},
                timeout=15,
            )
            if resp.status_code == 200:
                valid.append(key)
            else:
                console.log(
                    f"⚠️  Dropping invalid NCBI API key (...{suffix}): HTTP {resp.status_code}"
                )
        except Exception as e:
            console.log(f"⚠️  Dropping NCBI API key (...{suffix}): {e}")
    return valid


def download_contig_lengths(
    api_key: str | None = None,
    workers: int = 10,
    skip_assembly_summary: bool = False,
    api_keys: list[str] | None = None,
    per_key_concurrency: int = 3,
    requests_per_second: float | None = None,
    page_size: int = 100,
    resume: bool = True,
):
    """Fetch missing contig-length/sequence metadata via the official NCBI
    Datasets v2 ``sequence_reports`` API (see
    :mod:`hoodini.download.ncbi_sequence_reports`).
    """
    global NCBI_API_KEY
    NCBI_API_KEY = api_key
    if api_keys:
        api_keys = _validate_api_keys(api_keys)
        if not api_keys:
            console.log(
                "⚠️  No valid NCBI API keys remained after validation; "
                "falling back to single-session mode."
            )
        else:
            console.log(
                f"🔑 Using {len(api_keys)} valid NCBI API key(s), "
                f"{per_key_concurrency} connections each "
                f"(~{len(api_keys) * per_key_concurrency} concurrent requests)."
            )
    if not skip_assembly_summary:
        console.log("🔄 Updating local assembly_summary.parquet...")
        download_assembly_summary_db()
    else:
        console.log("⏭️  Skipping assembly_summary refresh (using local copy)")

    missing_df, latest_mtime = get_missing_contigs_from_summary(ASSEMBLY_SUMMARY)

    # Date-based filtering using remote file's Last-Modified date
    if missing_df.height > 0:
        # Get the remote file's Last-Modified date
        remote_last_mod = _get_remote_parquet_last_modified()

        if remote_last_mod:
            try:
                # Use DuckDB to check schema and filter by date
                con = duckdb.connect(":memory:")
                con.execute('SET memory_limit = "4GB"')

                # Check if seq_rel_date column exists
                schema_result = con.execute(
                    f"""
                    SELECT name FROM parquet_schema('{str(ASSEMBLY_SUMMARY)}')
                    WHERE name = 'seq_rel_date'
                """
                ).fetchone()

                if schema_result:
                    remote_date = remote_last_mod.date()
                    remote_date_str = remote_date.isoformat()

                    con.register("missing_asm", missing_df)

                    # Filter: keep if date is null OR date > remote_date
                    missing_df = con.execute(
                        f"""
                        SELECT DISTINCT m.assembly_accession
                        FROM missing_asm m
                        LEFT JOIN (
                            SELECT 
                                CAST(assembly_accession AS VARCHAR) as assembly_accession,
                                TRY_CAST(seq_rel_date AS DATE) as seq_rel_date
                            FROM read_parquet('{str(ASSEMBLY_SUMMARY)}')
                        ) a ON m.assembly_accession = a.assembly_accession
                        WHERE a.seq_rel_date IS NULL OR a.seq_rel_date > '{remote_date_str}'
                    """
                    ).pl()

                    console.log(
                        f"Filtered by remote date ({remote_date}): {missing_df.height:,} assemblies remain"
                    )
                else:
                    console.log(
                        "⚠️  'seq_rel_date' not found in assembly_summary; skipping date filtering"
                    )

                con.close()

            except Exception as e:
                console.log(f"⚠️  DuckDB date filtering failed: {e}; skipping date filter")

    if missing_df.height == 0:
        console.log("✅ No missing contig lengths to download.")
        return

    # Convert to list only at the end when needed for API call
    missing_list = missing_df["assembly_accession"].to_list()

    stats = asyncio.run(
        fetch_sequence_reports(
            missing_list,
            output_dir=CONTIG_LENGTHS_DIR,
            checkpoint_path=CONTIG_LENGTHS_CHECKPOINT,
            api_keys=api_keys,
            concurrency=workers,
            connections_per_key=per_key_concurrency,
            requests_per_second=requests_per_second,
            page_size=page_size,
            retries=MAX_RETRIES,
            target_mb=30,
            resume=resume,
        )
    )

    if stats.sequences_fetched == 0:
        console.log("✅ No contig length records returned.")
        return

    console.log(
        f"✅ Fetched sequence reports for {stats.assemblies_completed}/{len(missing_list)} "
        f"assemblies (skipped={stats.assemblies_skipped}, failed={stats.assemblies_failed})."
    )
    console.log(
        f"    sequences written={stats.sequences_fetched}; requests={stats.requests_made}; "
        f"retries={stats.retries}; HTTP 429s={stats.http_429}"
    )
    if stats.errors:
        console.log(f"⚠️  {len(stats.errors)} assemblies failed, e.g.: {stats.errors[:5]}")
