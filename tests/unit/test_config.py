"""
Unit tests for configuration loading and validation.

Tests:
- defaults.toml loading
- RuntimeConfig validation
- Config merging (defaults + TOML + CLI)
"""

from hoodini.config import load_default_config
from hoodini.config.settings import RuntimeConfig


class TestDefaultConfig:
    """Tests for default configuration loading."""

    def test_load_default_config_returns_dict(self):
        """Default config should return a dictionary."""
        defaults = load_default_config()
        assert isinstance(defaults, dict)

    def test_default_config_has_required_keys(self):
        """Default config should contain essential keys."""
        defaults = load_default_config()
        # Check some essential keys exist
        expected_keys = ["wn", "mod", "tree_mode", "cand_mode", "clust_method"]
        for key in expected_keys:
            assert key in defaults, f"Missing required key: {key}"

    def test_default_window_size(self):
        """Default window size should be reasonable."""
        defaults = load_default_config()
        assert defaults.get("wn", 0) > 0
        assert defaults.get("wn", 0) <= 100000  # Sanity check


class TestRuntimeConfig:
    """Tests for RuntimeConfig dataclass."""

    def test_runtime_config_minimal(self, tmp_path):
        """RuntimeConfig should work with minimal required fields."""
        cfg = RuntimeConfig(
            input_path="test.txt",
            output=str(tmp_path / "output"),
        )
        assert cfg.input_path == "test.txt"
        assert "output" in cfg.output

    def test_runtime_config_defaults(self, tmp_path):
        """RuntimeConfig should have sensible defaults."""
        cfg = RuntimeConfig(
            input_path="test.txt",
            output=str(tmp_path / "output"),
        )
        assert cfg.force is False
        assert cfg.debug is False
        assert cfg.keep is False
        # num_threads may be None (uses system default)

    def test_runtime_config_tree_modes(self, tmp_path):
        """RuntimeConfig should accept valid tree modes."""
        valid_modes = [
            "taxonomy",
            "aai_tree",
            "ani_tree",
            "fast_ml",
            "fast_nj",
            "neigh_phylo_tree",
            "neigh_similarity_tree",
            "foldmason_tree",
        ]
        for mode in valid_modes:
            cfg = RuntimeConfig(
                input_path="test.txt",
                output=str(tmp_path / "output"),
                tree_mode=mode,
            )
            assert cfg.tree_mode == mode

    def test_runtime_config_cand_modes(self, tmp_path):
        """RuntimeConfig should accept valid candidate modes."""
        valid_modes = ["any_ipg", "best_ipg", "best_id", "one_id", "same_id"]
        for mode in valid_modes:
            cfg = RuntimeConfig(
                input_path="test.txt",
                output=str(tmp_path / "output"),
                cand_mode=mode,
            )
            assert cfg.cand_mode == mode

    def test_runtime_config_window_modes(self, tmp_path):
        """RuntimeConfig should accept valid window modes."""
        for mod in ["win_nts", "win_genes"]:
            cfg = RuntimeConfig(
                input_path="test.txt",
                output=str(tmp_path / "output"),
                mod=mod,
            )
            assert cfg.mod == mod

    def test_runtime_config_annotations(self, tmp_path):
        """RuntimeConfig should accept annotation flags."""
        cfg = RuntimeConfig(
            input_path="test.txt",
            output=str(tmp_path / "output"),
            padloc=True,
            deffinder=True,
            cctyper=True,
            genomad=True,
            emapper=True,
            sorfs=True,
        )
        assert cfg.padloc is True
        assert cfg.deffinder is True
        assert cfg.cctyper is True
        assert cfg.genomad is True
        assert cfg.emapper is True
        assert cfg.sorfs is True


class TestConfigMerging:
    """Tests for configuration merging logic."""

    def test_cli_overrides_defaults(self, tmp_path):
        """CLI values should override defaults."""
        defaults = load_default_config()
        cli_values = {"wn": 25000, "tree_mode": "aai_tree"}

        merged = {**defaults, **cli_values}

        assert merged["wn"] == 25000
        assert merged["tree_mode"] == "aai_tree"

    def test_none_values_dont_override(self, tmp_path):
        """None values should not override defaults."""
        defaults = load_default_config()
        cli_values = {"wn": None, "tree_mode": "aai_tree"}

        # Filter out None values before merging
        cli_filtered = {k: v for k, v in cli_values.items() if v is not None}
        merged = {**defaults, **cli_filtered}

        assert merged["wn"] == defaults["wn"]
        assert merged["tree_mode"] == "aai_tree"
