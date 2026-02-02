"""
Pytest configuration and shared fixtures for Hoodini tests.
"""

import os
import sys
from pathlib import Path

import pytest

# Ensure hoodini is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Test data paths
TESTS_DIR = Path(__file__).parent
FEATURE_INPUTS = TESTS_DIR / "feature_coverage" / "inputs"


@pytest.fixture
def sample_protein_ids():
    """Sample protein IDs for testing."""
    return [
        "WP_217844005.1",
        "WP_347132630.1",
        "WP_239738697.1",
    ]


@pytest.fixture
def sample_assembly_ids():
    """Sample assembly IDs for testing."""
    return [
        "GCF_000005845.2",  # E. coli K-12
        "GCF_000009045.1",  # Bacillus subtilis
    ]


@pytest.fixture
def sample_nucleotide_ids():
    """Sample nucleotide accessions for testing."""
    return [
        "NC_000913.3",
        "NZ_CP028116.1",
    ]


@pytest.fixture
def tmp_output(tmp_path):
    """Temporary output directory for tests."""
    output = tmp_path / "output"
    output.mkdir()
    return output
