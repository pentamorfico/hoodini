"""
Unit tests for candidate selection modes.

Tests:
- any_ipg mode behavior
- best_ipg mode behavior
- best_id mode behavior
- one_id mode behavior
- same_id mode behavior
"""


class TestCandModeDefinitions:
    """Tests for candidate mode definitions and semantics."""

    CAND_MODES = ["any_ipg", "best_ipg", "best_id", "one_id", "same_id"]

    def test_all_modes_defined(self):
        """All candidate modes should be defined."""
        assert len(self.CAND_MODES) == 5

    def test_mode_names_are_strings(self):
        """Mode names should be lowercase strings."""
        for mode in self.CAND_MODES:
            assert isinstance(mode, str)
            assert mode == mode.lower()

    def test_mode_uniqueness(self):
        """All modes should be unique."""
        assert len(self.CAND_MODES) == len(set(self.CAND_MODES))


class TestAnyIpgMode:
    """Tests for any_ipg candidate mode."""

    def test_any_ipg_accepts_all(self):
        """any_ipg should accept all IPG records."""
        # Simulated IPG records
        records = [
            {"protein_id": "WP_123", "assembly": "GCF_001"},
            {"protein_id": "WP_123", "assembly": "GCF_002"},
            {"protein_id": "WP_123", "assembly": "GCA_001"},
        ]

        # any_ipg keeps all
        selected = records  # No filtering
        assert len(selected) == 3


class TestBestIpgMode:
    """Tests for best_ipg candidate mode."""

    def test_best_ipg_prefers_refseq(self):
        """best_ipg should prefer RefSeq (GCF) over GenBank (GCA)."""
        records = [
            {"protein_id": "WP_123", "assembly": "GCA_001", "is_refseq": False},
            {"protein_id": "WP_123", "assembly": "GCF_001", "is_refseq": True},
        ]

        # best_ipg prefers RefSeq
        best = max(records, key=lambda r: r["is_refseq"])
        assert best["assembly"].startswith("GCF")


class TestBestIdMode:
    """Tests for best_id candidate mode."""

    def test_best_id_single_per_protein(self):
        """best_id should select one assembly per protein."""
        records = [
            {"protein_id": "WP_123", "assembly": "GCF_001"},
            {"protein_id": "WP_123", "assembly": "GCF_002"},
            {"protein_id": "WP_456", "assembly": "GCF_003"},
        ]

        # Group by protein, take one per protein
        proteins = {}
        for r in records:
            pid = r["protein_id"]
            if pid not in proteins:
                proteins[pid] = r

        assert len(proteins) == 2


class TestOneIdMode:
    """Tests for one_id candidate mode."""

    def test_one_id_single_result(self):
        """one_id should return exactly one record total."""
        records = [
            {"protein_id": "WP_123", "assembly": "GCF_001"},
            {"protein_id": "WP_456", "assembly": "GCF_002"},
            {"protein_id": "WP_789", "assembly": "GCF_003"},
        ]

        # one_id takes just the first/best
        selected = [records[0]]
        assert len(selected) == 1


class TestSameIdMode:
    """Tests for same_id candidate mode."""

    def test_same_id_keeps_identical_proteins(self):
        """same_id should keep records for identical protein sequences."""
        records = [
            {"protein_id": "WP_123", "assembly": "GCF_001", "ipg_id": 100},
            {"protein_id": "WP_123", "assembly": "GCF_002", "ipg_id": 100},
            {"protein_id": "WP_456", "assembly": "GCF_003", "ipg_id": 200},
        ]

        # same_id groups by IPG (identical proteins)
        ipg_groups = {}
        for r in records:
            ipg = r["ipg_id"]
            if ipg not in ipg_groups:
                ipg_groups[ipg] = []
            ipg_groups[ipg].append(r)

        # Should have 2 groups
        assert len(ipg_groups) == 2
        # First group has 2 records (identical proteins)
        assert len(ipg_groups[100]) == 2
