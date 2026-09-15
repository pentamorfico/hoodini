"""Fetch NCBI sequence-level metadata via the official Datasets v2 REST API.

This replaces the previous approach of calling an undocumented, reverse
engineered internal RPC (``prefetch_links.make_seqrep_url`` /
``fetch_to_queue`` in :mod:`hoodini.download.contig_lengths`) with the
documented, versioned endpoint:

    POST https://api.ncbi.nlm.nih.gov/datasets/v2/genome/sequence_reports
    body: {"accession": "<assembly accession>", "page_size": <int>, "page_token": "<opt>"}

Empirically confirmed response shape (no OpenAPI spec is publicly reachable,
so this was verified with live requests against a range of real assemblies):

* A successful response is ``{"reports": [...], "total_count": <int>,
  "next_page_token": "<opt>"}``. ``total_count`` is only present on the
  first page and can *undercount* the true total for large/fragmented
  genomes, so callers must paginate until ``next_page_token`` is absent and
  must never use ``total_count`` to decide when to stop.
* An unknown/invalid/malformed accession returns **HTTP 200** with body
  ``{}`` (no ``reports`` key) rather than a 4xx. Code must treat a missing
  ``reports`` key as "no results", not rely on the HTTP status.
* Each report record is one *replicon/scaffold*-level sequence (chromosome,
  plasmid, unplaced/unlocalized scaffold, patch, ...), not a raw sequencing
  contig. NCBI's ``contig_count`` in ``assembly_summary`` can be far larger
  than the number of ``sequence_reports`` records for a given assembly
  (e.g. GCF_000002875.2: contig_count=278, sequence_reports records=76) —
  the API aggregates raw contigs into scaffolds before reporting them. If a
  caller needs true sub-scaffold contig boundaries, this endpoint is not
  sufficient; the AGP file or full FASTA would be required instead. This
  is the exact same granularity the previous undocumented endpoint already
  returned, so this migration does not change the semantics of Hoodini's
  existing ``contig_lengths`` data — only the transport used to fetch it.
* There is no field named ``molecule_type``. The closest real field is
  ``assigned_molecule_location_type`` (e.g. ``"Chromosome"``, ``"Plasmid"``,
  ``"Mitochondrion"``).
* Real fields observed across many assemblies: ``assembly_accession``,
  ``assembly_unit``, ``assigned_molecule_location_type``, ``chr_name``,
  ``gc_count`` (string-typed integer), ``gc_percent``, ``genbank_accession``,
  ``length``, ``refseq_accession``, ``role``, ``sequence_name``,
  ``sort_order``, ``ucsc_style_name``, ``unlocalized_count``.
"""

from __future__ import annotations

import asyncio
import contextlib
import sqlite3
import time
from collections.abc import AsyncIterator, Iterable
from dataclasses import dataclass, field
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any

import aiohttp  # type: ignore[import]
import pyarrow as pa
import pyarrow.parquet as pq  # type: ignore[import]
from rich.progress import (  # type: ignore[import]
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from hoodini.utils.logging_utils import console

SEQUENCE_REPORTS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/sequence_reports"

# NCBI-observed rate limits: 5 req/s anonymous, 10 req/s with a valid API
# key, ~3 req/s with an invalid one. Used only as safe defaults; callers can
# override via ``requests_per_second``.
DEFAULT_RPS_ANONYMOUS = 4.0
DEFAULT_RPS_WITH_KEY = 8.0

RETRYABLE_STATUSES = {429, 500, 502, 503, 504}
# A genome could legitimately need many pages (page_size=100 -> ~8 pages for
# the human genome); this is only a defensive circuit breaker against an
# API/logic bug that would otherwise loop forever on a repeating token.
MAX_PAGES_PER_ASSEMBLY = 20_000


class RetryableHTTPError(RuntimeError):
    """Raised when a request exhausted its retry budget."""


class TokenBucketRateLimiter:
    """Real requests-per-second limiter (distinct from bounding concurrency).

    A ``asyncio.Semaphore`` only bounds how many requests are *in flight* at
    once; it says nothing about how many are *issued* per second. This is a
    classic token bucket: tokens regenerate continuously at ``rate`` per
    second up to ``capacity``, and ``acquire()`` waits until a token is
    available before letting the caller proceed.
    """

    def __init__(self, rate: float, capacity: float | None = None):
        if rate <= 0:
            raise ValueError("rate must be > 0")
        self.rate = rate
        self.capacity = capacity if capacity is not None else max(rate, 1.0)
        self._tokens = self.capacity
        self._last = time.monotonic()
        self._lock = asyncio.Lock()

    async def acquire(self) -> None:
        async with self._lock:
            while True:
                now = time.monotonic()
                elapsed = now - self._last
                self._last = now
                self._tokens = min(self.capacity, self._tokens + elapsed * self.rate)
                if self._tokens >= 1:
                    self._tokens -= 1
                    return
                wait = (1 - self._tokens) / self.rate
                await asyncio.sleep(wait)


def _parse_retry_after(value: str | None) -> float | None:
    """Parse a ``Retry-After`` header: either delta-seconds or an HTTP-date."""
    if not value:
        return None
    value = value.strip()
    with contextlib.suppress(ValueError):
        return max(0.0, float(value))
    with contextlib.suppress(Exception):
        dt = parsedate_to_datetime(value)
        return max(0.0, (dt - dt.now(dt.tzinfo)).total_seconds())
    return None


@dataclass
class FetchStats:
    """Counters surfaced to callers/tests and used to drive progress output."""

    assemblies_completed: int = 0
    assemblies_failed: int = 0
    assemblies_skipped: int = 0
    sequences_fetched: int = 0
    requests_made: int = 0
    retries: int = 0
    http_429: int = 0
    errors: list[str] = field(default_factory=list)


class SequenceReportCheckpoint:
    """SQLite-backed record of assemblies that have been fully downloaded.

    Unlike the pre-existing anti-join against ``contig_lengths/*.parquet``
    (which can treat an assembly as "done" even if only some of its rows
    made it into a committed part file before an interruption), an
    assembly is only marked done here once *every* page of its
    sequence_reports has been fetched and its rows handed to the writer.
    """

    def __init__(self, path: Path):
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(path))
        self._conn.execute(
            "CREATE TABLE IF NOT EXISTS completed_assemblies ("
            "assembly_accession TEXT PRIMARY KEY, "
            "sequences_written INTEGER NOT NULL, "
            "completed_at TEXT NOT NULL)"
        )
        self._conn.execute(
            "CREATE TABLE IF NOT EXISTS failed_assemblies ("
            "assembly_accession TEXT PRIMARY KEY, "
            "error TEXT, "
            "failed_at TEXT NOT NULL)"
        )
        self._conn.commit()

    def is_done(self, accession: str) -> bool:
        row = self._conn.execute(
            "SELECT 1 FROM completed_assemblies WHERE assembly_accession = ?",
            (accession,),
        ).fetchone()
        return row is not None

    def mark_done(self, accession: str, sequences_written: int) -> None:
        self._conn.execute(
            "INSERT OR REPLACE INTO completed_assemblies VALUES (?, ?, datetime('now'))",
            (accession, sequences_written),
        )
        self._conn.execute(
            "DELETE FROM failed_assemblies WHERE assembly_accession = ?", (accession,)
        )
        self._conn.commit()

    def mark_failed(self, accession: str, error: str) -> None:
        self._conn.execute(
            "INSERT OR REPLACE INTO failed_assemblies VALUES (?, ?, datetime('now'))",
            (accession, error[:500]),
        )
        self._conn.commit()

    def completed_count(self) -> int:
        (n,) = self._conn.execute("SELECT COUNT(*) FROM completed_assemblies").fetchone()
        return int(n)

    def close(self) -> None:
        self._conn.close()


class PartRotatingWriter:
    """Write buffered rows to rotating ``part-NNNNN.parquet`` files.

    Column set is the union of every row seen so far (rows are dicts with
    heterogeneous optional keys); the four backward-compatible columns
    (``genbankAccession``, ``refseqAccession``, ``assemblyAccession``,
    ``length``) always exist so downstream readers relying on that schema
    keep working even for partitions that lack, say, RefSeq accessions.
    """

    CORE_TYPES: dict[str, pa.DataType] = {
        "genbankAccession": pa.string(),
        "refseqAccession": pa.string(),
        "assemblyAccession": pa.string(),
        "length": pa.int64(),
    }

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
        columns = dict.fromkeys(self.CORE_TYPES)
        for row in rows:
            columns.update(dict.fromkeys(row))
        table = pa.table(
            {
                name: pa.array([row.get(name) for row in rows], type=self.CORE_TYPES.get(name))
                for name in columns
            }
        )
        tmp = self.dir / f"part-{self.part_idx:05d}.parquet.tmp"
        pq.write_table(table, tmp, compression="zstd")
        return tmp.stat().st_size

    def _commit_tmp(self) -> None:
        tmp = self.dir / f"part-{self.part_idx:05d}.parquet.tmp"
        final = self.dir / f"part-{self.part_idx:05d}.parquet"
        tmp.replace(final)
        self.part_idx += 1
        self.total_files += 1

    def add_many(self, rows: list[dict[str, Any]]) -> None:
        if not rows:
            return
        self._buffer.extend(rows)
        self.total_rows += len(rows)
        while len(self._buffer) >= self.rows_target:
            self._flush_targeted()

    def _flush_targeted(self) -> None:
        if not self._buffer:
            return
        rows = self._buffer
        target_n = min(self.rows_target, len(rows))
        try_rows = rows[:target_n]
        self._write_once(try_rows)
        self._commit_tmp()
        del rows[:target_n]
        self._buffer = rows

    def close(self) -> None:
        if self._buffer:
            self._write_once(self._buffer)
            self._commit_tmp()
            self._buffer.clear()


def normalize_sequence_report(raw: dict[str, Any], accession: str) -> dict[str, Any]:
    """Map a real v2 ``sequence_reports`` record onto Hoodini's row schema.

    The first four keys are the pre-existing backward-compatible contract
    (:data:`hoodini.utils.contig_metadata.CONTIG_SCHEMA`); everything else is
    new, additive metadata. ``sequenceId`` is a stable per-row identifier:
    RefSeq accession if present, else GenBank accession, else a
    ``assembly:sequence_name`` fallback for the rare record with neither.
    """
    genbank = raw.get("genbank_accession")
    refseq = raw.get("refseq_accession")
    sequence_name = raw.get("sequence_name")
    sequence_id = refseq or genbank or (f"{accession}:{sequence_name}" if sequence_name else None)

    gc_count_raw = raw.get("gc_count")
    gc_count: int | None
    try:
        gc_count = int(gc_count_raw) if gc_count_raw is not None else None
    except (TypeError, ValueError):
        gc_count = None

    return {
        "genbankAccession": genbank,
        "refseqAccession": refseq,
        "assemblyAccession": raw.get("assembly_accession") or accession,
        "length": raw.get("length"),
        "sequenceId": sequence_id,
        "sequenceName": sequence_name,
        "chrName": raw.get("chr_name"),
        "role": raw.get("role"),
        "assemblyUnit": raw.get("assembly_unit"),
        "assignedMoleculeLocationType": raw.get("assigned_molecule_location_type"),
        "gcPercent": raw.get("gc_percent"),
        "gcCount": gc_count,
        "sortOrder": raw.get("sort_order"),
    }


async def _post_with_retries(
    session: aiohttp.ClientSession,
    body: dict[str, Any],
    rate_limiter: TokenBucketRateLimiter,
    stats: FetchStats,
    retries: int,
    timeout: float,
) -> dict[str, Any]:
    attempt = 0
    while True:
        await rate_limiter.acquire()
        stats.requests_made += 1
        try:
            async with session.post(
                SEQUENCE_REPORTS_URL, json=body, timeout=aiohttp.ClientTimeout(total=timeout)
            ) as resp:
                if resp.status == 429:
                    stats.http_429 += 1
                if resp.status in RETRYABLE_STATUSES:
                    if attempt >= retries:
                        raise RetryableHTTPError(
                            f"HTTP {resp.status} for {body.get('accession')} after {attempt} retries"
                        )
                    delay = _parse_retry_after(resp.headers.get("Retry-After"))
                    if delay is None:
                        delay = min(2**attempt, 30)
                    stats.retries += 1
                    attempt += 1
                    await asyncio.sleep(delay)
                    continue
                resp.raise_for_status()
                return await resp.json()
        except (
            TimeoutError,
            aiohttp.ClientConnectionError,
            aiohttp.ServerDisconnectedError,
        ):
            if attempt >= retries:
                raise
            stats.retries += 1
            attempt += 1
            await asyncio.sleep(min(2**attempt, 30))


async def fetch_sequence_report_page(
    session: aiohttp.ClientSession,
    accession: str,
    page_size: int,
    page_token: str | None,
    rate_limiter: TokenBucketRateLimiter,
    stats: FetchStats,
    retries: int = 5,
    timeout: float = 30.0,
) -> dict[str, Any]:
    """Fetch one page of ``sequence_reports`` for ``accession``."""
    body: dict[str, Any] = {"accession": accession, "page_size": page_size}
    if page_token:
        body["page_token"] = page_token
    return await _post_with_retries(session, body, rate_limiter, stats, retries, timeout)


async def fetch_assembly_sequences(
    session: aiohttp.ClientSession,
    accession: str,
    page_size: int,
    rate_limiter: TokenBucketRateLimiter,
    stats: FetchStats,
    retries: int = 5,
    timeout: float = 30.0,
) -> AsyncIterator[list[dict[str, Any]]]:
    """Yield each page's ``reports`` list for ``accession``, in order.

    Pagination relies solely on ``next_page_token`` presence: ``total_count``
    is only returned on the first page and has been observed to undercount
    the true total for large/fragmented genomes.
    """
    page_token: str | None = None
    for _ in range(MAX_PAGES_PER_ASSEMBLY):
        data = await fetch_sequence_report_page(
            session, accession, page_size, page_token, rate_limiter, stats, retries, timeout
        )
        reports = data.get("reports") or []
        if reports:
            yield reports
        page_token = data.get("next_page_token")
        if not page_token:
            return


async def _worker(
    accession_queue: asyncio.Queue[str | None],
    rows_queue: asyncio.Queue[list[dict[str, Any]] | None],
    session: aiohttp.ClientSession,
    rate_limiter: TokenBucketRateLimiter,
    checkpoint: SequenceReportCheckpoint | None,
    stats: FetchStats,
    page_size: int,
    retries: int,
    timeout: float,
    resume: bool,
    progress_cb: Any,
) -> None:
    while True:
        accession = await accession_queue.get()
        if accession is None:
            return
        if resume and checkpoint is not None and checkpoint.is_done(accession):
            stats.assemblies_skipped += 1
            progress_cb()
            continue
        row_count = 0
        try:
            async for page_rows in fetch_assembly_sequences(
                session, accession, page_size, rate_limiter, stats, retries, timeout
            ):
                normalized = [normalize_sequence_report(r, accession) for r in page_rows]
                row_count += len(normalized)
                await rows_queue.put(normalized)
            stats.sequences_fetched += row_count
            stats.assemblies_completed += 1
            if checkpoint is not None:
                checkpoint.mark_done(accession, row_count)
        except Exception as exc:  # noqa: BLE001 - one bad assembly must not kill the run
            stats.assemblies_failed += 1
            stats.errors.append(f"{accession}: {exc}")
            if checkpoint is not None:
                checkpoint.mark_failed(accession, str(exc))
        finally:
            progress_cb()


async def _rows_consumer(
    rows_queue: asyncio.Queue[list[dict[str, Any]] | None],
    n_workers: int,
    writer: PartRotatingWriter,
) -> None:
    done = 0
    while done < n_workers:
        item = await rows_queue.get()
        if item is None:
            done += 1
            continue
        writer.add_many(item)
    writer.close()


async def fetch_sequence_reports(
    accessions: Iterable[str],
    output_dir: Path,
    *,
    checkpoint_path: Path | None = None,
    api_keys: list[str] | None = None,
    concurrency: int = 20,
    connections_per_key: int = 3,
    requests_per_second: float | None = None,
    page_size: int = 100,
    retries: int = 5,
    timeout: float = 30.0,
    target_mb: float = 30.0,
    resume: bool = True,
    show_progress: bool = True,
) -> FetchStats:
    """Fetch sequence-level metadata for many assemblies concurrently.

    With ``api_keys`` given, one aiohttp session per key is created (each
    sending an ``api-key`` header and its own rate-limit bucket), with
    ``connections_per_key`` workers per session — mirroring the multi-key
    fan-out already used for the legacy endpoint. Without keys, a single
    anonymous session runs ``concurrency`` workers. ``requests_per_second``
    is applied *per session* (i.e. per key, since NCBI's rate limit is a
    per-key/per-IP bucket); it defaults to a conservative value for
    anonymous vs. keyed access based on observed NCBI limits (5 req/s
    anonymous, 10 req/s keyed).
    """
    stats = FetchStats()
    writer = PartRotatingWriter(output_dir, target_bytes=int(target_mb * 1024 * 1024))
    checkpoint = SequenceReportCheckpoint(checkpoint_path) if checkpoint_path else None

    accession_list = list(dict.fromkeys(accessions))
    if not accession_list:
        writer.close()
        if checkpoint is not None:
            checkpoint.close()
        return stats

    keys: list[str | None] = list(api_keys) if api_keys else [None]
    n_workers_per_session = connections_per_key if keys != [None] else concurrency
    total_workers = n_workers_per_session * len(keys)
    rps = requests_per_second or (DEFAULT_RPS_WITH_KEY if keys != [None] else DEFAULT_RPS_ANONYMOUS)

    accession_queue: asyncio.Queue[str | None] = asyncio.Queue()
    for acc in accession_list:
        accession_queue.put_nowait(acc)
    for _ in range(total_workers):
        accession_queue.put_nowait(None)

    rows_queue: asyncio.Queue[list[dict[str, Any]] | None] = asyncio.Queue(maxsize=50)
    consumer_task = asyncio.create_task(_rows_consumer(rows_queue, total_workers, writer))

    progress = Progress(
        SpinnerColumn(),
        BarColumn(),
        TextColumn("{task.completed}/{task.total} assemblies"),
        TimeElapsedColumn(),
        console=console,
        disable=not show_progress,
    )
    task_id = progress.add_task("Fetching sequence_reports", total=len(accession_list))

    def progress_cb() -> None:
        progress.update(task_id, advance=1)

    try:
        async with contextlib.AsyncExitStack() as stack:
            worker_tasks = []
            for key in keys:
                headers = {"Accept": "application/json", "User-Agent": "hoodini-sequence-reports"}
                if key:
                    headers["api-key"] = key
                connector = aiohttp.TCPConnector(limit_per_host=64, ttl_dns_cache=300)
                session = await stack.enter_async_context(
                    aiohttp.ClientSession(headers=headers, connector=connector)
                )
                rate_limiter = TokenBucketRateLimiter(rate=rps)
                for _ in range(n_workers_per_session):
                    worker_tasks.append(
                        asyncio.create_task(
                            _worker(
                                accession_queue,
                                rows_queue,
                                session,
                                rate_limiter,
                                checkpoint,
                                stats,
                                page_size,
                                retries,
                                timeout,
                                resume,
                                progress_cb,
                            )
                        )
                    )
            with progress:
                await asyncio.gather(*worker_tasks)
            for _ in range(total_workers):
                await rows_queue.put(None)
            await consumer_task
    finally:
        if checkpoint is not None:
            checkpoint.close()

    return stats
