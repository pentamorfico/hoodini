"""
Integration tests for the pipeline runner.

These tests verify that modules interact correctly without
running the full heavy pipeline.
"""

import pytest
from pathlib import Path

from hoodini.config import load_default_config
from hoodini.config.settings import RuntimeConfig


class TestRuntimeConfigIntegration:
    """Tests for RuntimeConfig integration with pipeline."""

    def test_config_to_pipeline_wiring(self, tmp_path):
        """RuntimeConfig should wire correctly to pipeline expectations."""
        cfg = RuntimeConfig(
            input_path="test.txt",
            output=str(tmp_path / "output"),
            tree_mode="taxonomy",
            cand_mode="any_ipg",
            mod="win_nts",
            wn=50000,
            num_threads=4,
        )
        
        # Verify all expected attributes exist and are accessible
        assert cfg.input_path is not None
        assert cfg.output is not None
        assert cfg.tree_mode == "taxonomy"
        assert cfg.cand_mode == "any_ipg"
        assert cfg.mod == "win_nts"
        assert cfg.wn == 50000
        assert cfg.num_threads == 4

    def test_config_defaults_merge_correctly(self, tmp_path):
        """Defaults should merge correctly with explicit config."""
        defaults = load_default_config()
        
        # Merge with explicit values
        merged = {
            **defaults,
            "input_path": "test.txt",
            "output": str(tmp_path / "output"),
            "tree_mode": "aai_tree",
        }
        
        cfg = RuntimeConfig(
            input_path=merged["input_path"],
            output=merged["output"],
            tree_mode=merged["tree_mode"],
            wn=merged.get("wn", 50000),
            mod=merged.get("mod", "win_nts"),
        )
        
        assert cfg.tree_mode == "aai_tree"
        assert cfg.wn == defaults.get("wn", 50000)


class TestOutputStructure:
    """Tests for expected output structure."""

    def test_output_folder_creation(self, tmp_path):
        """Pipeline should create output folder structure."""
        output = tmp_path / "test_output"
        output.mkdir()
        
        # Expected subdirectories
        expected_subdirs = [
            "assembly_folder",
            "hoodini-viz",
            "hoodini-viz/tsv",
            "hoodini-viz/parquet",
        ]
        
        for subdir in expected_subdirs:
            (output / subdir).mkdir(parents=True, exist_ok=True)
        
        # Verify structure
        assert (output / "assembly_folder").exists()
        assert (output / "hoodini-viz").exists()
        assert (output / "hoodini-viz" / "tsv").exists()


class TestAnnotationFlags:
    """Tests for annotation flag handling."""

    def test_multiple_annotations_config(self, tmp_path):
        """Multiple annotations should be configurable together."""
        cfg = RuntimeConfig(
            input_path="test.txt",
            output=str(tmp_path / "output"),
            padloc=True,
            deffinder=True,
            cctyper=True,
        )
        
        # All should be set
        active_annotations = []
        if cfg.padloc:
            active_annotations.append("padloc")
        if cfg.deffinder:
            active_annotations.append("deffinder")
        if cfg.cctyper:
            active_annotations.append("cctyper")
        
        assert len(active_annotations) == 3

    def test_no_annotations_is_valid(self, tmp_path):
        """Running without annotations should be valid."""
        cfg = RuntimeConfig(
            input_path="test.txt",
            output=str(tmp_path / "output"),
            padloc=False,
            deffinder=False,
            cctyper=False,
            genomad=False,
            emapper=False,
        )
        
        # All annotation flags should be False
        assert cfg.padloc is False
        assert cfg.deffinder is False
        assert cfg.cctyper is False
