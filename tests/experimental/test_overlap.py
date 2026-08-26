"""Experimental checks on the fork-BWA benchmark BAM.

The benchmark fixture lives under tests/benchmark/ and is not shipped, so these
tests skip cleanly when it is absent (rather than breaking test collection).
"""
import os

import pytest
import pysam

_BAM = "tests/benchmark/verify_final_perf.bam"


@pytest.fixture(scope="module")
def bam():
    if not os.path.exists(_BAM):
        pytest.skip("benchmark BAM not present")
    return pysam.AlignmentFile(_BAM, "rb")


def test_fork_overlap_counts(bam):
    """Count the fork-BWA origin (ST) tags; only checks it parses cleanly."""
    ori1_primary = 0
    ori2_secondary = 0
    for read in bam:
        if not read.has_tag("ST"):
            continue
        st = read.get_tag("ST")
        if not read.is_secondary and st == 1:
            ori1_primary += 1
        elif read.is_secondary and st == 2:
            ori2_secondary += 1
    assert ori1_primary >= 0 and ori2_secondary >= 0
