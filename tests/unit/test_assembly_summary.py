"""Regression coverage for typed and atomic assembly-summary updates."""

from pathlib import Path

import polars as pl
import pytest

from hoodini.download import assembly_summary

HEADER = [
    "assembly_accession",
    "taxid",
    "species_taxid",
    "organism_name",
    "ftp_path",
    "group",
    "isolate",
    "genome_size",
    "gc_percent",
]


def report(accession="GCA_000000001.1", *, header=None, rows=None, prefix="#"):
    header = HEADER if header is None else header
    default = {
        "assembly_accession": accession,
        "taxid": "2",
        "species_taxid": "2",
        "organism_name": "Test bacterium",
        "ftp_path": "https://example.test/assembly",
        "group": "bacteria",
        "isolate": "001",
        "genome_size": "1000",
        "gc_percent": "42.5",
    }
    rows = [default] if rows is None else [{**default, **row} for row in rows]
    return (
        "# Assembly report fixture\n"
        + prefix
        + "\t".join(header)
        + "\n"
        + "".join("\t".join(str(row.get(col, "na")) for col in header) + "\n" for row in rows)
    )


@pytest.fixture
def serve_reports(monkeypatch):
    calls = []

    def configure(contents):
        def download(urls, dest_dir, **kwargs):
            directory = Path(dest_dir)
            directory.mkdir(parents=True, exist_ok=True)
            calls.append(directory)
            paths = []
            for name, content in zip(kwargs["out_names"], contents):
                if content is None:
                    continue
                path = directory / name
                path.write_text(content)
                paths.append(str(path))
            return paths

        monkeypatch.setattr(assembly_summary, "download_with_aria2c", download)
        return calls

    return configure


def update(output, *, historical=False, columns=None):
    assembly_summary.download_assembly_db(
        ["genbank"], output, columns_to_keep=columns, include_historical=historical
    )


@pytest.mark.parametrize("column", ["group", "isolate", "pubmed_id"])
def test_text_columns_do_not_depend_on_first_thousand_rows(tmp_path, serve_reports, column):
    header = [*HEADER, "pubmed_id"]
    rows = [{"assembly_accession": f"GCA_{i:09d}.1", column: "001"} for i in range(1, 1101)]
    rows.append({"assembly_accession": "GCA_000001101.1", column: "bacteria"})
    serve_reports([report(header=header, rows=rows)])
    output = tmp_path / "summary.parquet"
    update(output)
    frame = pl.read_parquet(output)
    assert frame.height == 1101
    assert frame.schema[column] == pl.Utf8
    assert set(frame[column]) == {"001", "bacteria"}
    assert frame.schema["taxid"] == pl.Int64
    assert frame.schema["genome_size"] == pl.Int64
    assert frame.schema["gc_percent"] == pl.Float64


@pytest.mark.parametrize("prefix", ["#", "# "])
def test_header_names_define_types_after_column_reordering(tmp_path, serve_reports, prefix):
    header = ["new_column", *reversed(HEADER)]
    serve_reports([report(header=header, prefix=prefix, rows=[{"new_column": "001"}])])
    output = tmp_path / "summary.parquet"
    update(output)
    row = pl.read_parquet(output).row(0, named=True)
    assert row["assembly_accession"] == "GCA_000000001.1"
    assert row["group"] == "bacteria"
    assert row["genome_size"] == 1000
    assert row["new_column"] == "001"


def test_current_and_historical_reports_can_have_different_optional_columns(
    tmp_path, serve_reports
):
    serve_reports(
        [
            report(header=[*HEADER, "future_field"], rows=[{"future_field": "001"}]),
            report("GCA_000000002.1", header=[name for name in HEADER if name != "gc_percent"]),
        ]
    )
    output = tmp_path / "summary.parquet"
    update(output, historical=True)
    frame = pl.read_parquet(output).sort("assembly_accession")
    assert frame.height == 2
    assert frame["future_field"].to_list() == ["001", None]
    assert frame["gc_percent"].to_list() == [42.5, None]


def test_requested_optional_fields_remain_present_with_typed_nulls(tmp_path, serve_reports):
    serve_reports([report(header=[name for name in HEADER if name != "gc_percent"])])
    output = tmp_path / "summary.parquet"
    update(output, columns=["assembly_accession", "gc_percent", "future_field"])
    frame = pl.read_parquet(output)
    assert frame.columns == ["assembly_accession", "gc_percent", "future_field"]
    assert frame["gc_percent"].to_list() == [None]
    assert frame.schema["gc_percent"] == pl.Float64
    assert frame.schema["future_field"] == pl.Utf8


@pytest.mark.parametrize(
    "bad_report",
    [
        report("GCA_000000002.1", rows=[{"genome_size": "bacteria"}]),
        report("GCA_000000002.1").rsplit("\t", 1)[0] + "\n",
        report("GCA_000000002.1").rstrip("\n") + "\textra\n",
        "#assembly_accession\ttaxid\n",
        "not an assembly summary\n",
        report(rows=[{"assembly_accession": "na"}]),
        report(header=[*HEADER, "taxid"]),
    ],
    ids=[
        "wrong_numeric_type",
        "short_row",
        "extra_field",
        "empty",
        "no_header",
        "no_id",
        "duplicate_header",
    ],
)
def test_bad_source_cannot_replace_existing_database(tmp_path, serve_reports, bad_report):
    calls = serve_reports([report(), bad_report])
    output = tmp_path / "summary.parquet"
    pl.DataFrame({"assembly_accession": ["GCA_999999999.1"]}).write_parquet(output)
    original = output.read_bytes()
    with pytest.raises((ValueError, RuntimeError)):
        update(output, historical=True)
    assert output.read_bytes() == original
    assert (calls[0] / "assembly_summary_genbank_historical.txt").read_text() == bad_report


@pytest.mark.parametrize("contents", [[report(), None], []])
def test_incomplete_download_cannot_publish_a_partial_database(tmp_path, serve_reports, contents):
    serve_reports(contents)
    output = tmp_path / "summary.parquet"
    output.write_bytes(b"previous database")
    with pytest.raises((FileNotFoundError, RuntimeError)):
        update(output, historical=True)
    assert output.read_bytes() == b"previous database"


def test_failed_parquet_write_preserves_previous_database(tmp_path, serve_reports, monkeypatch):
    calls = serve_reports([report()])
    output = tmp_path / "summary.parquet"
    output.write_bytes(b"previous database")

    def failed_write(self, path, **kwargs):
        Path(path).write_bytes(b"incomplete parquet")
        raise OSError("simulated disk full")

    monkeypatch.setattr(pl.DataFrame, "write_parquet", failed_write)
    with pytest.raises(OSError, match="simulated disk full"):
        update(output)
    assert output.read_bytes() == b"previous database"
    assert calls[0] != tmp_path
    assert (calls[0] / "assembly_summary_genbank.txt").is_file()


def test_success_publishes_all_sources_and_cleans_staging(tmp_path, serve_reports):
    calls = serve_reports([report(f"GCA_{i:09d}.1") for i in range(1, 5)])
    stale_file = tmp_path / "assembly_summary_refseq.txt"
    stale_file.write_text("stale download")
    output = tmp_path / "summary.parquet"
    assembly_summary.download_assembly_db(["refseq", "genbank"], output)
    assert pl.read_parquet(output).height == 4
    assert stale_file.read_text() == "stale download"
    assert calls[0] != tmp_path
    assert not calls[0].exists()


def test_invalid_first_install_does_not_create_output(tmp_path, serve_reports):
    serve_reports([report(rows=[{"genome_size": "bacteria"}])])
    output = tmp_path / "summary.parquet"
    with pytest.raises(ValueError):
        update(output)
    assert not output.exists()


def test_initialization_stops_before_overwriting_results_when_summary_setup_fails(
    tmp_path, monkeypatch
):
    from hoodini.pipeline import initialize

    root = tmp_path / "package"
    (root / "data").mkdir(parents=True)
    results = tmp_path / "results"
    results.mkdir()
    marker = results / "previous-result.txt"
    marker.write_text("keep these results")

    def invalid_summary(output_path):
        raise ValueError("invalid assembly summary")

    def must_not_continue():
        pytest.fail("Pipeline continued after assembly-summary setup failed")

    monkeypatch.setattr(initialize, "files", lambda package: root)
    monkeypatch.setattr(assembly_summary, "download_assembly_summary_db", invalid_summary)
    monkeypatch.setattr(initialize, "check_contig_lengths_db", must_not_continue)
    with pytest.raises(ValueError, match="invalid assembly summary"):
        initialize.initialize_inputs(input_path="query.txt", output=results, force=True)
    assert marker.read_text() == "keep these results"
