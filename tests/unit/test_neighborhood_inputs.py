"""Regressions for local neighborhood selection and mixed input formats."""

from pathlib import Path

import polars as pl
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from hoodini.pipeline.helpers.neighborhood_extractor import extract_neighborhood


@pytest.fixture
def local_assembly(tmp_path):
    gff = tmp_path / "assembly.gff"
    gff.write_text(
        "ctg1\tlocal\tCDS\t10\t39\t2.5\t+\t0\tID=p1\n"
        "ctg2\tlocal\tCDS\t10\t39\t3.5\t+\t0\tID=p2\n"
        "ctg1\tlocal\tCDS\t50\t79\t4.5\t+\t0\tID=p3\n"
        "ctg2\tlocal\tCDS\t50\t79\t5.5\t+\t0\tID=p4\n"
    )
    faa = tmp_path / "assembly.faa"
    faa.write_text("".join(f">p{i}\n{'M' * 10}\n" for i in range(1, 5)))
    fna = tmp_path / "assembly.fna"
    fna.write_text(f">ctg1\n{'A' * 100}\n>ctg2\n{'C' * 150}\n")
    return {"gff_file": str(gff), "faa_file": str(faa), "fna_file": str(fna)}


def extract_local(local_assembly, **kwargs):
    params = {
        "gbf_file": None,
        "protein_id": None,
        "nucleotide_id": "ctg1",
        "input_type": "nucleotide",
        "unique_id": "hood1",
        "window": 100,
        **local_assembly,
        **kwargs,
    }
    return extract_neighborhood(**params)


@pytest.mark.parametrize("with_fna", [True, False])
@pytest.mark.parametrize("mode", ["win_nts", "win_genes"])
def test_nucleotide_neighborhood_excludes_other_contigs(local_assembly, with_fna, mode):
    if not with_fna:
        local_assembly["fna_file"] = None
    genes, hood, _, error = extract_local(local_assembly, mode=mode, start=10, end=79)
    assert error is None
    assert genes["protein_id"].to_list() == ["p1", "p3"]
    assert genes["seqid"].unique().to_list() == ["ctg1"]
    assert hood["seqid"].to_list() == ["ctg1"]


def test_protein_query_infers_contig_and_matching_fna(local_assembly):
    genes, hood, _, error = extract_local(
        local_assembly, nucleotide_id=None, protein_id="p2", input_type="protein"
    )
    assert error is None
    assert genes["protein_id"].to_list() == ["p2", "p4"]
    assert hood["seqid"].to_list() == ["ctg2"]
    assert set(hood["sequence"][0]) == {"C"}


@pytest.mark.parametrize("nucleotide_id", ["missing", None])
def test_missing_or_ambiguous_contig_reports_selection_error(local_assembly, nucleotide_id):
    local_assembly["fna_file"] = None
    genes, hood, _, error = extract_local(local_assembly, nucleotide_id=nucleotide_id)
    assert genes is None
    assert hood is None
    assert error
    assert "contig" in error.lower() or "nucleotide" in error.lower()


def test_mixed_gbff_and_numeric_gff_can_be_concatenated(local_assembly, tmp_path):
    record = SeqRecord(Seq("A" * 100), id="ctg1.1", name="ctg1", description="fixture")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(9, 39, strand=1),
            type="CDS",
            qualifiers={"protein_id": ["gb1"], "translation": ["M" * 10]},
        )
    ]
    gbff = tmp_path / "assembly.gbff"
    SeqIO.write(record, gbff, "genbank")
    local_genes, _, _, local_error = extract_local(local_assembly)
    gb_genes, _, _, gb_error = extract_local(
        local_assembly, gbf_file=str(gbff), nucleotide_id="ctg1.1"
    )
    assert local_error is None
    assert gb_error is None
    combined = pl.concat([gb_genes, local_genes], how="diagonal")
    assert combined.height == 3
    assert combined.schema["score"] == pl.Utf8
    assert combined.schema["phase"] == pl.Utf8


def test_gff_comments_are_skipped_before_parsing(local_assembly):
    gff = Path(local_assembly["gff_file"])
    gff.write_text("##gff-version 3\n# short comment\n" + gff.read_text())
    genes, _, _, error = extract_local(local_assembly)
    assert error is None
    assert genes["protein_id"].to_list() == ["p1", "p3"]


def test_single_contig_input_still_works_without_an_accession(local_assembly):
    gff = Path(local_assembly["gff_file"])
    gff.write_text(
        "\n".join(line for line in gff.read_text().splitlines() if line.startswith("ctg2"))
    )
    genes, hood, _, error = extract_local(local_assembly, nucleotide_id=None)
    assert error is None
    assert genes["protein_id"].to_list() == ["p2", "p4"]
    assert hood["seqid"].to_list() == ["ctg2"]
    assert set(hood["sequence"][0]) == {"C"}


def test_protein_id_on_multiple_contigs_requires_nucleotide_id(local_assembly):
    gff = Path(local_assembly["gff_file"])
    gff.write_text(gff.read_text().replace("ID=p2", "ID=p1"))
    genes, hood, _, error = extract_local(
        local_assembly, nucleotide_id=None, protein_id="p1", input_type="protein"
    )
    assert genes is None
    assert hood is None
    assert "unique GFF contig" in error
