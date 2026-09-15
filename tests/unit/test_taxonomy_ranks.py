"""Exercise actual metadata assembly using deterministic NCBI lineage fixtures."""

import polars as pl
import pytest

from hoodini.pipeline import taxonomy


@pytest.mark.parametrize(
    ("ranks", "names", "expected"),
    [
        ({2: "domain"}, {2: "Bacteria"}, "Bacteria"),
        ({2: "superkingdom"}, {2: "Bacteria"}, "Bacteria"),
        (
            {2: "domain", 3: "superkingdom"},
            {2: "modern", 3: "legacy"},
            "legacy",
        ),
        (
            {2: "superkingdom", 3: "domain"},
            {2: "legacy", 3: "modern"},
            "legacy",
        ),
        ({2: "domain"}, {2: "Archaea"}, "Archaea"),
        ({2: "acellular root"}, {2: "Viruses"}, "unclassified"),
        ({}, {}, "unclassified"),
    ],
)
def test_top_rank_compatibility(monkeypatch, ranks, names, expected):
    class TaxonomyFixture:
        def get_lineage_translator(self, taxids):
            return {tid: [1, *ranks] for tid in taxids}

        def get_rank(self, taxids):
            return {1: "cellular root", **ranks}

        def get_taxid_translator(self, taxids):
            return {1: "cellular organisms", **names}

    monkeypatch.setattr(taxonomy, "NCBITaxa", TaxonomyFixture)
    records = pl.DataFrame({"taxid": [2], "unique_id": ["1"], "og_index": [0], "failed": [None]})
    neighborhoods = pl.DataFrame(
        {
            "unique_id": ["1"],
            "start_win": [0],
            "end_win": [100],
            "strand_win": ["+"],
            "start_target": [10],
            "end_target": [40],
        }
    )
    result = taxonomy._build_leaf_metadata(records, neighborhoods)
    assert result["superkingdom"].to_list() == [expected]
    assert result["species"].to_list() == ["unclassified"]
    assert "domain" not in result.columns
