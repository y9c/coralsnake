"""Shared fixtures for the coralsnake test suite."""

from pathlib import Path

import pytest


@pytest.fixture
def data_dir():
    """Path to tests/data/."""
    return Path(__file__).resolve().parent / "data"


@pytest.fixture
def has_integration_data(data_dir):
    """Skip if core integration data (ref + reads) is missing."""
    required = ["ref.fa", "test1.fq", "test2.fq"]
    if not all((data_dir / f).exists() for f in required):
        pytest.skip("integration test data not present")


@pytest.fixture
def has_liftover_data(data_dir):
    """Skip if liftover test data is missing."""
    required = ["liftover_input.bam", "liftover_annot.tsv", "liftover_faidx.fai"]
    if not all((data_dir / f).exists() for f in required):
        pytest.skip("liftover test data not present")


@pytest.fixture
def has_gtf_data(data_dir):
    """Skip if GTF/FASTA test data is missing."""
    required = ["R64-1-1.release57.gtf", "R64-1-1.fa"]
    if not all((data_dir / f).exists() for f in required):
        pytest.skip("GTF test data not present")
