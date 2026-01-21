import asyncio
import contextlib
import json
import os
from datetime import datetime
from email.utils import parsedate_to_datetime
from importlib.resources import files
from pathlib import Path
from typing import Any

import aiohttp  # type: ignore[import]
import duckdb
import polars as pl  # type: ignore[import]
import pyarrow as pa
import pyarrow.parquet as pq  # type: ignore[import]
import requests  # type: ignore[import]
from rich.progress import (  # type: ignore[import]
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from hoodini.download.assembly_summary import download_assembly_summary_db
from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table
from hoodini.utils.logging_utils import console

NCBI_API_KEY: str | None = os.environ.get("NCBI_API_KEY")

DATA_DIR = files("hoodini").joinpath("data")
CONTIG_LENGTHS_DIR = DATA_DIR.joinpath("contig_lengths")
MASTER_CONTIGS = DATA_DIR.joinpath("contig_lengths.parquet")
ASSEMBLY_SUMMARY = DATA_DIR.joinpath("assembly_summary.parquet")

DEFAULT_GROUPS = {"bacteria", "viral", "archaea", "metagenomes", "other"}
MAX_RETRIES = 3
REMOTE_CONTIG_LENGTHS_URL = "https://storage.hoodini.bio/contig_lengths.parquet"

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
        con = duckdb.connect(":memory:")
        con.execute('SET memory_limit = "4GB"')

        # Create temp table for allowed assemblies if provided
        if allowed_assemblies_df is not None:
            allowed_list = allowed_assemblies_df["assembly_accession"].to_list()
            con.execute("CREATE TEMP TABLE allowed_asm (assembly_accession VARCHAR)")
            con.executemany("INSERT INTO allowed_asm VALUES (?)", [(a,) for a in allowed_list])

        # Build the query for valid assemblies from assembly_summary
        groups_str = ", ".join(f"'{g}'" for g in DEFAULT_GROUPS)
        summary_query = f"""
            SELECT DISTINCT CAST(assembly_accession AS VARCHAR) as assembly_accession
            FROM read_parquet('{str(ASSEMBLY_SUMMARY)}')
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
            missing_df = con.execute(
                f"""
                WITH summary AS ({summary_query}),
                existing AS (
                    SELECT DISTINCT CAST(assemblyAccession AS VARCHAR) as assembly_accession
                    FROM read_parquet('{contig_glob}')
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

        con.close()

    except Exception as e:
        console.log(f"⚠️  DuckDB failed, falling back to Polars streaming: {e}")
        # Fallback to Polars if DuckDB fails
        summary_lf = (
            pl.scan_parquet(str(ASSEMBLY_SUMMARY))
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
                pl.scan_parquet(
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


SENTINEL = object()


def _consume_lines(buffer: bytearray):
    while True:
        nl = buffer.find(b"\n")
        if nl == -1:
            break
        line = buffer[:nl]
        del buffer[: nl + 1]
        yield line.decode("utf-8", errors="ignore").strip()


class PartRotatingWriter:
    def __init__(
        self, dataset_dir: Path, target_bytes: int = 30 * 1024 * 1024, start_rows: int = 80_000
    ):
        self.dir = dataset_dir
        self.dir.mkdir(parents=True, exist_ok=True)
        self.target_bytes = target_bytes
        self.rows_target = start_rows
        self._buffer: list[dict[str, Any]] = []
        self.part_idx = self._next_index()
        self.total_rows = 0
        self.total_files = 0

    def _next_index(self) -> int:
        existing = list(self.dir.glob("part-*.parquet"))
        if not existing:
            return 0
        max_id = -1
        for p in existing:
            with contextlib.suppress(Exception):
                max_id = max(max_id, int(p.stem.split("-")[1]))
        return max_id + 1

    def _write_once(self, rows: list[dict[str, Any]]) -> int:
        table = pa.Table.from_pylist(rows)
        tmp = self.dir / f"part-{self.part_idx:05d}.parquet.tmp"
        pq.write_table(table, tmp, compression="zstd")
        return tmp.stat().st_size

    def _commit_tmp(self):
        tmp = self.dir / f"part-{self.part_idx:05d}.parquet.tmp"
        final = self.dir / f"part-{self.part_idx:05d}.parquet"
        tmp.replace(final)
        self.part_idx += 1
        self.total_files += 1

    def add_many(self, rows: list[dict[str, Any]]):
        if not rows:
            return
        self._buffer.extend(rows)
        self.total_rows += len(rows)
        while len(self._buffer) >= self.rows_target:
            self._flush_targeted()

    def _flush_targeted(self):
        if not self._buffer:
            return
        rows = self._buffer
        target_n = min(self.rows_target, len(rows))
        try_rows = rows[:target_n]
        self._write_once(try_rows)
        self._commit_tmp()
        del rows[:target_n]
        self._buffer = rows

    def close(self):
        if self._buffer:
            _ = self._write_once(self._buffer)
            self._commit_tmp()
            self._buffer.clear()


async def fetch_to_queue(session, asm, url, queue, batch_rows, retries, timeout_s, sem):
    attempt = 0
    CHUNK = 1 << 19
    sent = 0
    async with sem:
        while True:
            attempt += 1
            try:
                batch: list[dict[str, Any]] = []
                async with session.get(url, timeout=timeout_s) as resp:
                    if resp.status != 200:
                        raise RuntimeError(f"HTTP {resp.status}")
                    buf = bytearray()
                    async for chunk in resp.content.iter_chunked(CHUNK):
                        if not chunk:
                            continue
                        buf.extend(chunk)
                        for line in _consume_lines(buf):
                            if not line:
                                continue
                            try:
                                obj = json.loads(line)
                            except Exception:
                                continue
                            batch.append(obj)
                            if len(batch) >= batch_rows:
                                await queue.put(batch)
                                sent += len(batch)
                                batch = []
                    if buf:
                        tail = buf.decode("utf-8", errors="ignore").strip()
                        if tail:
                            try:
                                obj = json.loads(tail)
                                batch.append(obj)
                            except Exception:
                                pass
                if batch:
                    await queue.put(batch)
                    sent += len(batch)
                await queue.put(SENTINEL)
                return True, f"rows={sent}", sent
            except Exception as e:
                if attempt <= retries:
                    await asyncio.sleep(min(2**attempt, 10))
                    continue
                await queue.put(SENTINEL)
                return False, str(e), sent


async def writer_consumer(queue, n_producers, writer):
    done = 0
    while True:
        item = await queue.get()
        if item is SENTINEL:
            done += 1
            if done >= n_producers:
                break
            continue
        writer.add_many(item)
    writer.close()
    return writer.total_rows, writer.total_files


async def stream_and_write(
    pairs, target_mb=30, batch_rows=5000, concurrency=10, retries=3, timeout=60
):
    target_bytes = target_mb * 1024 * 1024
    writer = PartRotatingWriter(dataset_dir=CONTIG_LENGTHS_DIR, target_bytes=target_bytes)
    queue = asyncio.Queue(maxsize=20)
    sem = asyncio.Semaphore(concurrency)
    consumer_task = asyncio.create_task(
        writer_consumer(queue, n_producers=len(pairs), writer=writer)
    )

    timeout_cfg = aiohttp.ClientTimeout(total=None, connect=30, sock_read=600)
    headers = {"Accept-Encoding": "gzip, deflate", "User-Agent": "seqrep-dl/hoodini"}
    connector = aiohttp.TCPConnector(limit_per_host=64, ttl_dns_cache=300)

    async with aiohttp.ClientSession(
        timeout=timeout_cfg, connector=connector, headers=headers
    ) as session:
        progress = Progress(
            SpinnerColumn(),
            BarColumn(),
            TextColumn("{task.completed}/{task.total} {task.description}"),
            TimeElapsedColumn(),
        )
        task_id = progress.add_task(f"Downloading (conc={concurrency})", total=len(pairs))
        ok = 0
        total_rows = 0
        with progress:
            tasks = [
                asyncio.create_task(
                    fetch_to_queue(session, asm, url, queue, batch_rows, retries, timeout, sem)
                )
                for asm, url in pairs
            ]
            for fut in asyncio.as_completed(tasks):
                ok1, info, sent = await fut
                total_rows += sent
                progress.update(task_id, advance=1, description=info[:120])
                if ok1:
                    ok += 1
    total_rows_written, files_written = await consumer_task
    return ok, total_rows, total_rows_written, files_written


def download_contig_lengths(
    api_key: str | None = None, workers: int = 10, skip_assembly_summary: bool = False
):
    global NCBI_API_KEY
    NCBI_API_KEY = api_key
    if not skip_assembly_summary:
        console.log("🔄 Updating local assembly_summary.parquet...")
        download_assembly_summary_db()
    else:
        console.log("⏭️  Skipping assembly_summary refresh (using local copy)")

    # Build allowed_assemblies as a DataFrame (not a Python set)
    allowed_assemblies_df: pl.DataFrame | None = None
    try:
        # Use DuckDB to get candidate IDs without loading full DataFrame
        con = duckdb.connect(":memory:")
        con.execute('SET memory_limit = "4GB"')
        groups_str = ", ".join(f"'{g}'" for g in DEFAULT_GROUPS)
        candidate_df = con.execute(
            f"""
            SELECT DISTINCT CAST(assembly_accession AS VARCHAR) as assembly_accession
            FROM read_parquet('{str(ASSEMBLY_SUMMARY)}')
            WHERE "group" IN ({groups_str})
              AND ftp_path IS NOT NULL
              AND TRIM(ftp_path) != ''
              AND LOWER(ftp_path) != 'na'
        """
        ).pl()
        con.close()
        candidate_ids = candidate_df["assembly_accession"].to_list()

        if candidate_ids:
            links_df = get_prefetched_link_table(candidate_ids, kinds=["sequence_report"])
            # Keep as DataFrame instead of converting to Python set
            allowed_assemblies_df = (
                links_df.filter(pl.col("filetype") == "sequence_report")
                .select(pl.col("assembly_id").cast(pl.Utf8).alias("assembly_accession"))
                .unique()
            )
    except Exception:
        allowed_assemblies_df = None

    missing_df, latest_mtime = get_missing_contigs_from_summary(
        ASSEMBLY_SUMMARY, allowed_assemblies_df=allowed_assemblies_df
    )

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
                    SELECT column_name FROM parquet_schema('{str(ASSEMBLY_SUMMARY)}')
                    WHERE column_name = 'seq_rel_date'
                """
                ).fetchone()

                if schema_result:
                    remote_date = remote_last_mod.date()
                    remote_date_str = remote_date.isoformat()

                    # Create temp table for missing assemblies
                    missing_list = missing_df["assembly_accession"].to_list()
                    con.execute("CREATE TEMP TABLE missing_asm (assembly_accession VARCHAR)")
                    con.executemany(
                        "INSERT INTO missing_asm VALUES (?)", [(a,) for a in missing_list]
                    )

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
    df_links = get_prefetched_link_table(missing_list, kinds=["sequence_report"], seqrep_only=True)

    # Use Polars filter instead of boolean mask indexing
    links_filtered = df_links.filter(pl.col("filetype") == "sequence_report")
    pairs: list[tuple[str, str]] = list(
        zip(links_filtered["assembly_id"].to_list(), links_filtered["url"].to_list())
    )

    if not pairs:
        console.log("✅ No sequence_report links available for missing assemblies.")
        return

    ok, rows_fetched, rows_written, files_written = asyncio.run(
        stream_and_write(
            pairs,
            target_mb=30,
            batch_rows=5000,
            concurrency=workers,
            retries=MAX_RETRIES,
            timeout=60,
        )
    )

    if rows_written == 0:
        console.log("✅ No contig length records returned.")
        return

    console.log(
        f"✅ Downloaded contig lengths for {len(pairs)} sequence_report links (ok={ok}/{len(pairs)})."
    )
    console.log(
        f"    rows fetched={rows_fetched}; rows written={rows_written}; new parts={files_written}"
    )
