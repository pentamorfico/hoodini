"""Tests for the official NCBI Datasets v2 sequence_reports fetcher.

Uses lightweight fakes for aiohttp (no new test dependency): a session whose
``.post()`` returns a queued response or raises a queued exception, used as
``async with session.post(...) as resp``, matching real aiohttp usage.
"""

import asyncio

import pyarrow.parquet as pq
import pytest

from hoodini.download import ncbi_sequence_reports as nsr


class FakeResponse:
    def __init__(self, status=200, json_data=None, headers=None):
        self.status = status
        self._json = json_data if json_data is not None else {}
        self.headers = headers or {}

    async def __aenter__(self):
        return self

    async def __aexit__(self, *exc):
        return False

    async def json(self):
        return self._json

    def raise_for_status(self):
        if self.status >= 400:
            raise RuntimeError(f"HTTP {self.status}")


class FakeSession:
    """Async-context-manager session backed by a shared per-accession queue."""

    def __init__(self, responses: dict[str, list]):
        self._responses = responses
        self.posts: list[dict] = []

    async def __aenter__(self):
        return self

    async def __aexit__(self, *exc):
        return False

    def post(self, url, json, timeout=None):
        self.posts.append(json)
        acc = json["accession"]
        queue = self._responses[acc]
        item = queue.pop(0)
        if isinstance(item, BaseException):
            raise item
        return item


class FakeConnector:
    def __init__(self, *a, **k):
        pass

    async def close(self):
        pass


def _patch_aiohttp(monkeypatch, responses: dict[str, list]):
    """Route every ``aiohttp.ClientSession(...)`` call to a FakeSession sharing
    ``responses``, and make ``aiohttp.TCPConnector`` a no-op."""
    monkeypatch.setattr(nsr.aiohttp, "ClientSession", lambda **kw: FakeSession(responses))
    monkeypatch.setattr(nsr.aiohttp, "TCPConnector", FakeConnector)


@pytest.fixture(autouse=True)
def _no_real_sleep(monkeypatch):
    """Skip real backoff delays; we only need to exercise the retry logic."""

    async def fast_sleep(_delay):
        return None

    monkeypatch.setattr(nsr.asyncio, "sleep", fast_sleep)


def report(**overrides):
    base = {
        "assembly_accession": "GCF_000005845.2",
        "genbank_accession": "CP009273.1",
        "refseq_accession": "NC_000913.3",
        "length": 4641652,
        "gc_percent": 50.79,
        "gc_count": "2357528",
        "role": "assembled-molecule",
        "sequence_name": "ANW",
        "chr_name": "",
        "assembly_unit": "Primary Assembly",
    }
    base.update(overrides)
    return base


def make_rate_limiter():
    return nsr.TokenBucketRateLimiter(rate=1000.0)


# --- unit-level tests -------------------------------------------------


def test_normalize_sequence_report_maps_real_fields():
    row = nsr.normalize_sequence_report(report(), "GCF_000005845.2")
    assert row["assemblyAccession"] == "GCF_000005845.2"
    assert row["genbankAccession"] == "CP009273.1"
    assert row["refseqAccession"] == "NC_000913.3"
    assert row["length"] == 4641652
    assert row["sequenceId"] == "NC_000913.3"
    assert row["gcCount"] == 2357528
    assert row["role"] == "assembled-molecule"


def test_normalize_sequence_report_falls_back_when_no_accessions():
    row = nsr.normalize_sequence_report({"length": 10, "sequence_name": "scaffold1"}, "GCF_111.1")
    assert row["sequenceId"] == "GCF_111.1:scaffold1"
    assert row["assemblyAccession"] == "GCF_111.1"


@pytest.mark.parametrize(
    "value,expected",
    [
        (None, None),
        ("5", 5.0),
        ("Wed, 21 Oct 2015 07:28:00 GMT", None),  # only checked for type, not exact value
    ],
)
def test_parse_retry_after_delta_seconds(value, expected):
    result = nsr._parse_retry_after(value)
    if expected is not None:
        assert result == expected
    elif value is None:
        assert result is None
    else:
        assert isinstance(result, float)


def test_checkpoint_marks_done_and_resumes(tmp_path):
    cp = nsr.SequenceReportCheckpoint(tmp_path / "checkpoint.sqlite")
    assert not cp.is_done("GCF_1")
    cp.mark_done("GCF_1", 3)
    assert cp.is_done("GCF_1")
    assert cp.completed_count() == 1
    cp.close()

    # Reopening from disk preserves the completed set (resumability).
    cp2 = nsr.SequenceReportCheckpoint(tmp_path / "checkpoint.sqlite")
    assert cp2.is_done("GCF_1")
    cp2.close()


# --- fetch_assembly_sequences (pagination) -----------------------------


def test_single_page_assembly(monkeypatch):
    responses = {"GCF_A": [FakeResponse(200, {"reports": [report()], "total_count": 1})]}
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        pages = [
            page
            async for page in nsr.fetch_assembly_sequences(
                session, "GCF_A", 100, make_rate_limiter(), stats
            )
        ]
        return pages, stats

    pages, stats = asyncio.run(run())
    assert len(pages) == 1
    assert len(pages[0]) == 1
    assert stats.requests_made == 1


def test_multi_page_assembly(monkeypatch):
    responses = {
        "GCF_B": [
            FakeResponse(200, {"reports": [report()], "next_page_token": "tok1"}),
            FakeResponse(200, {"reports": [report(sequence_name="p2")], "next_page_token": "tok2"}),
            FakeResponse(200, {"reports": [report(sequence_name="p3")]}),
        ]
    }
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        pages = [
            page
            async for page in nsr.fetch_assembly_sequences(
                session, "GCF_B", 1, make_rate_limiter(), stats
            )
        ]
        return pages, stats

    pages, stats = asyncio.run(run())
    assert len(pages) == 3
    assert stats.requests_made == 3
    # total_count is intentionally never trusted for stopping pagination.
    assert sum(len(p) for p in pages) == 3


def test_429_then_retry_after_then_success(monkeypatch):
    responses = {
        "GCF_C": [
            FakeResponse(429, {}, headers={"Retry-After": "1"}),
            FakeResponse(200, {"reports": [report()]}),
        ]
    }
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        pages = [
            page
            async for page in nsr.fetch_assembly_sequences(
                session, "GCF_C", 100, make_rate_limiter(), stats, retries=3
            )
        ]
        return pages, stats

    pages, stats = asyncio.run(run())
    assert len(pages) == 1
    assert stats.http_429 == 1
    assert stats.retries == 1


def test_transient_500_then_success(monkeypatch):
    responses = {
        "GCF_D": [
            FakeResponse(500, {}),
            FakeResponse(200, {"reports": [report()]}),
        ]
    }
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        pages = [
            page
            async for page in nsr.fetch_assembly_sequences(
                session, "GCF_D", 100, make_rate_limiter(), stats, retries=3
            )
        ]
        return pages, stats

    pages, stats = asyncio.run(run())
    assert len(pages) == 1
    assert stats.retries == 1


def test_timeout_then_success(monkeypatch):
    responses = {
        "GCF_E": [
            TimeoutError("simulated timeout"),
            FakeResponse(200, {"reports": [report()]}),
        ]
    }
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        pages = [
            page
            async for page in nsr.fetch_assembly_sequences(
                session, "GCF_E", 100, make_rate_limiter(), stats, retries=3
            )
        ]
        return pages, stats

    pages, stats = asyncio.run(run())
    assert len(pages) == 1
    assert stats.retries == 1


def test_invalid_accession_returns_no_pages(monkeypatch):
    # NCBI returns HTTP 200 with an empty body (no "reports" key) for an
    # unknown/invalid accession -- never a 4xx.
    responses = {"GCF_BOGUS": [FakeResponse(200, {})]}
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        pages = [
            page
            async for page in nsr.fetch_assembly_sequences(
                session, "GCF_BOGUS", 100, make_rate_limiter(), stats
            )
        ]
        return pages

    pages = asyncio.run(run())
    assert pages == []


def test_empty_reports_list_stops_pagination(monkeypatch):
    responses = {"GCF_F": [FakeResponse(200, {"reports": []})]}
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        pages = [
            page
            async for page in nsr.fetch_assembly_sequences(
                session, "GCF_F", 100, make_rate_limiter(), stats
            )
        ]
        return pages

    assert asyncio.run(run()) == []


def test_exhausted_retries_raise(monkeypatch):
    responses = {"GCF_G": [FakeResponse(500, {}), FakeResponse(500, {})]}
    _patch_aiohttp(monkeypatch, responses)

    async def run():
        session = nsr.aiohttp.ClientSession()
        stats = nsr.FetchStats()
        async for _ in nsr.fetch_assembly_sequences(
            session, "GCF_G", 100, make_rate_limiter(), stats, retries=1
        ):
            pass

    with pytest.raises(nsr.RetryableHTTPError):
        asyncio.run(run())


# --- fetch_sequence_reports (end-to-end orchestration) -----------------


def test_fetch_sequence_reports_partial_failure_and_resume(tmp_path, monkeypatch):
    responses = {
        "GCF_OK1": [FakeResponse(200, {"reports": [report(assembly_accession="GCF_OK1")]})],
        "GCF_OK2": [FakeResponse(200, {"reports": [report(assembly_accession="GCF_OK2")]})],
        "GCF_BAD": [FakeResponse(500, {}), FakeResponse(500, {})],
    }
    _patch_aiohttp(monkeypatch, responses)

    out_dir = tmp_path / "out"
    checkpoint_path = tmp_path / "checkpoint.sqlite"

    stats = asyncio.run(
        nsr.fetch_sequence_reports(
            ["GCF_OK1", "GCF_OK2", "GCF_BAD"],
            output_dir=out_dir,
            checkpoint_path=checkpoint_path,
            concurrency=2,
            retries=1,
            show_progress=False,
        )
    )

    assert stats.assemblies_completed == 2
    assert stats.assemblies_failed == 1
    assert stats.sequences_fetched == 2
    assert len(stats.errors) == 1

    # Incremental Parquet write: rows for the two successful assemblies made
    # it to disk despite the third assembly's total failure.
    parts = sorted(out_dir.glob("part-*.parquet"))
    assert parts
    table = pq.read_table(parts[0])
    assert set(table["assemblyAccession"].to_pylist()) == {"GCF_OK1", "GCF_OK2"}

    # Resume: re-running with the same checkpoint skips the two completed
    # assemblies and only re-attempts the failed one.
    responses["GCF_BAD"] = [FakeResponse(200, {"reports": [report(assembly_accession="GCF_BAD")]})]
    stats2 = asyncio.run(
        nsr.fetch_sequence_reports(
            ["GCF_OK1", "GCF_OK2", "GCF_BAD"],
            output_dir=out_dir,
            checkpoint_path=checkpoint_path,
            concurrency=2,
            retries=1,
            show_progress=False,
        )
    )
    assert stats2.assemblies_skipped == 2
    assert stats2.assemblies_completed == 1
    assert stats2.assemblies_failed == 0


def test_fetch_sequence_reports_empty_input_is_a_noop(tmp_path):
    stats = asyncio.run(
        nsr.fetch_sequence_reports([], output_dir=tmp_path / "out", show_progress=False)
    )
    assert stats.assemblies_completed == 0
    assert stats.sequences_fetched == 0
