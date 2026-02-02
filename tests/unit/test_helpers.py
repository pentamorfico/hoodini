"""
Unit tests for helper utilities.

Tests:
- nuc2asmlen functionality
- prefetch_links generation
- Other utility functions
"""

import pytest


class TestNuc2Asmlen:
    """Tests for nuc2asmlen utility."""

    def test_nuc2asmlen_with_list(self, sample_nucleotide_ids):
        """nuc2asmlen should accept a list of accessions."""
        from hoodini.pipeline.helpers.nuc2asmlen import run_nuc2asmlen
        
        df = run_nuc2asmlen(sample_nucleotide_ids)
        
        assert df is not None
        assert hasattr(df, "height")  # Polars DataFrame
        # Should have expected columns
        assert "NucleotideAccession" in df.columns
        assert "AssemblyAccession" in df.columns
        assert "length" in df.columns

    def test_nuc2asmlen_returns_polars(self, sample_nucleotide_ids):
        """nuc2asmlen should return a Polars DataFrame."""
        import polars as pl
        from hoodini.pipeline.helpers.nuc2asmlen import run_nuc2asmlen
        
        df = run_nuc2asmlen(sample_nucleotide_ids)
        
        assert isinstance(df, pl.DataFrame)


class TestPrefetchLinks:
    """Tests for prefetch_links utility."""

    def test_prefetch_links_generates_urls(self, sample_assembly_ids):
        """prefetch_links should generate valid URLs."""
        from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table
        
        df = get_prefetched_link_table(sample_assembly_ids, kinds=["gbff"])
        
        assert df.height > 0
        assert "url" in df.columns
        
        # URLs should start with https
        urls = df["url"].to_list()
        for url in urls:
            assert url.startswith("https://")

    def test_prefetch_links_multiple_kinds(self, sample_assembly_ids):
        """prefetch_links should handle multiple file kinds."""
        from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table
        
        kinds = ["gbff", "gff", "fna"]
        df = get_prefetched_link_table(sample_assembly_ids, kinds=kinds)
        
        # Should have entries for each kind per assembly
        assert df.height == len(sample_assembly_ids) * len(kinds)

    def test_prefetch_links_sequence_report(self, sample_assembly_ids):
        """prefetch_links should handle sequence_report kind."""
        from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table
        
        df = get_prefetched_link_table(sample_assembly_ids, kinds=["sequence_report"])
        
        assert df.height == len(sample_assembly_ids)
        assert "sequence_report" in df["filetype"].to_list()


class TestAccessionParsing:
    """Tests for accession parsing utilities."""

    def test_refseq_assembly_pattern(self):
        """RefSeq assemblies should match GCF pattern."""
        import re
        
        refseq_ids = ["GCF_000005845.2", "GCF_000009045.1"]
        pattern = r"^GCF_\d+\.\d+$"
        
        for acc in refseq_ids:
            assert re.match(pattern, acc), f"{acc} should match RefSeq pattern"

    def test_genbank_assembly_pattern(self):
        """GenBank assemblies should match GCA pattern."""
        import re
        
        genbank_ids = ["GCA_001886335.1", "GCA_000001405.28"]
        pattern = r"^GCA_\d+\.\d+$"
        
        for acc in genbank_ids:
            assert re.match(pattern, acc), f"{acc} should match GenBank pattern"

    def test_protein_accession_pattern(self):
        """Protein accessions should match WP/NP/YP patterns."""
        import re
        
        protein_ids = ["WP_217844005.1", "NP_001234.1", "YP_005678.1"]
        pattern = r"^[WNY]P_\d+\.\d+$"
        
        for acc in protein_ids:
            assert re.match(pattern, acc), f"{acc} should match protein pattern"
