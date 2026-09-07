"""Regression checks for incremental NCBI contig metadata updates (#86)."""

import duckdb
import polars as pl
import pyarrow.parquet as pq
import pytest

from hoodini.download import contig_lengths
from hoodini.pipeline import parse_ipg
from hoodini.pipeline.helpers import nuc2asmlen


def force_polars_fallback(monkeypatch):
    def unavailable(*args, **kwargs):
        raise RuntimeError("Exercise the Polars fallback")

    monkeypatch.setattr(duckdb, "connect", unavailable)


@pytest.fixture
def metadata_dir(tmp_path, monkeypatch):
    root = tmp_path / "hoodini"
    data = root / "data"
    parts = data / "contig_lengths"
    parts.mkdir(parents=True)
    for module in (nuc2asmlen, parse_ipg):
        monkeypatch.setattr(module, "files", lambda package: root)
    monkeypatch.setattr(contig_lengths, "CONTIG_LENGTHS_DIR", parts)
    monkeypatch.setattr(contig_lengths, "ASSEMBLY_SUMMARY", data / "assembly_summary.parquet")
    return data, parts


def write_parts(parts, shape):
    rows = [
        {"genbankAccession": "GB1", "assemblyAccession": "GCA_1", "length": 100},
        {"genbankAccession": "GB2", "assemblyAccession": "GCF_2", "length": 200},
    ]
    if shape != "no_refseq":
        rows[1]["refseqAccession"] = "RS2"
    if shape == "refseq_first":
        rows.reverse()
    if shape == "refseq_only":
        rows = [{"refseqAccession": "RS2", "assemblyAccession": "GCF_2", "length": 200}]
    for i, row in enumerate(rows):
        pl.DataFrame([row]).write_parquet(parts / f"part-{i:05d}.parquet")


@pytest.mark.parametrize("fallback", [False, True])
@pytest.mark.parametrize("shape", ["refseq_first", "refseq_last", "no_refseq", "refseq_only"])
def test_lookup_handles_heterogeneous_or_absent_optional_columns(
    metadata_dir, monkeypatch, fallback, shape
):
    _, parts = metadata_dir
    write_parts(parts, shape)
    if fallback:
        force_polars_fallback(monkeypatch)
    result = nuc2asmlen.run_nuc2asmlen(["GB1", "RS2", "missing"])
    assert result["NucleotideAccession"].to_list() == ["GB1", "RS2", "missing"]
    assert result["length"].to_list() == [
        None if shape == "refseq_only" else 100,
        None if shape == "no_refseq" else 200,
        None,
    ]


@pytest.mark.parametrize("fallback", [False, True])
def test_missing_assembly_detection_tolerates_partition_schema_changes(
    metadata_dir, monkeypatch, fallback
):
    data, parts = metadata_dir
    write_parts(parts, "refseq_last")
    summary = data / "assembly_summary.parquet"
    pl.DataFrame(
        {
            "assembly_accession": ["GCA_1", "GCF_2", "GCA_3"],
            "group": ["bacteria"] * 3,
            "ftp_path": ["https://example.test/assembly"] * 3,
        }
    ).write_parquet(summary)
    if fallback:
        force_polars_fallback(monkeypatch)
    missing, mtime = contig_lengths.get_missing_contigs_from_summary(summary)
    assert missing["assembly_accession"].to_list() == ["GCA_3"]
    assert mtime is not None


def test_writer_keeps_optional_fields_from_later_rows_and_types(metadata_dir):
    _, parts = metadata_dir
    writer = contig_lengths.PartRotatingWriter(parts, start_rows=2)
    writer.add_many(
        [
            {"genbankAccession": "GB1", "assemblyAccession": "GCA_1", "length": 100},
            {
                "genbankAccession": "GB2",
                "refseqAccession": "RS2",
                "assemblyAccession": "GCF_2",
                "length": 200,
                "sequenceName": "second-row-only",
            },
            {"genbankAccession": "GB3", "assemblyAccession": "GCA_3", "length": 300},
        ]
    )
    writer.close()
    first = pq.read_table(parts / "part-00000.parquet")
    second = pq.read_table(parts / "part-00001.parquet")
    assert first["refseqAccession"].to_pylist() == [None, "RS2"]
    assert first["sequenceName"].to_pylist() == [None, "second-row-only"]
    assert second["refseqAccession"].to_pylist() == [None]
    assert first.schema.field("refseqAccession") == second.schema.field("refseqAccession")


def test_pipeline_contig_lookup_uses_normalized_scan(metadata_dir, monkeypatch):
    data, parts = metadata_dir
    write_parts(parts, "refseq_last")
    pl.DataFrame(
        {
            "assembly_accession": ["GCA_1"],
            "taxid": [2],
            "species_taxid": [2],
            "organism_name": ["Bacteria"],
            "infraspecific_name": [None],
            "assembly_level": ["Contig"],
            "group": ["bacteria"],
        }
    ).write_parquet(data / "assembly_summary.parquet")

    def no_fallback(accessions):
        pytest.fail("Direct lookup should succeed without falling back to nuc2asmlen")

    monkeypatch.setattr(parse_ipg, "run_nuc2asmlen", no_fallback)
    records = pl.DataFrame(
        {
            "nucleotide_id": ["GB1"],
            "input_type": ["nucleotide"],
            "failed": [None],
            "premade": [False],
            "start": [None],
            "end": [None],
        }
    )
    result = parse_ipg._fetch_nucleotide_data(records)
    assert result["assembly_id"].to_list() == ["GCA_1"]
    assert result["sequence_length"].to_list() == [100]
    assert result["taxid"].to_list() == [2]
