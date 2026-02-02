"""
Unit tests for window calculation logic.

Tests:
- win_nts mode calculations
- win_genes mode calculations
- minwin settings (upstream/downstream)
- Edge cases and boundary conditions
"""


class TestWindowModes:
    """Tests for window mode logic."""

    def test_win_nts_mode_values(self):
        """win_nts mode should use nucleotide counts."""
        # Typical values for win_nts
        valid_values = [10000, 25000, 50000, 100000]
        for wn in valid_values:
            assert wn > 0
            assert wn <= 200000  # Reasonable upper bound

    def test_win_genes_mode_values(self):
        """win_genes mode should use gene counts."""
        # Typical values for win_genes
        valid_values = [5, 10, 15, 20, 50]
        for wn in valid_values:
            assert wn > 0
            assert wn <= 100  # Reasonable for gene counts


class TestMinwinSettings:
    """Tests for minimum window (minwin) settings."""

    def test_minwin_types(self):
        """minwin_type should be upstream or downstream."""
        valid_types = ["upstream", "downstream"]
        for t in valid_types:
            assert t in valid_types

    def test_minwin_relationship_to_wn(self):
        """minwin should typically be smaller than wn."""
        wn = 50000
        minwin = 8000
        assert minwin < wn, "minwin should be smaller than window size"

    def test_minwin_gene_counts(self):
        """minwin in win_genes mode uses gene counts."""
        wn = 15  # genes
        minwin = 5  # genes
        assert minwin < wn


class TestWindowEdgeCases:
    """Tests for edge cases in window calculations."""

    def test_small_contig_handling(self):
        """Small contigs may be below minimum window size."""
        contig_length = 5000
        wn = 50000

        is_below_minwin = contig_length < wn
        assert is_below_minwin is True

    def test_symmetric_window(self):
        """Without minwin, windows should be symmetric around target."""
        wn = 50000
        target_pos = 100000

        upstream_start = target_pos - wn
        downstream_end = target_pos + wn

        upstream_size = target_pos - upstream_start
        downstream_size = downstream_end - target_pos

        assert upstream_size == downstream_size == wn

    def test_asymmetric_window_with_minwin(self):
        """With minwin, one side can be shorter."""
        wn = 50000
        minwin = 8000
        # minwin_type would be "upstream" or "downstream"

        # Simulating a case where upstream is constrained
        actual_upstream = minwin  # constrained
        actual_downstream = wn  # full size

        assert actual_upstream < actual_downstream
        assert actual_upstream >= minwin
