"""Experimental fork-BWA benchmark smoke test.

The benchmark index and ``bwamem`` may be absent, so these tests skip cleanly
(rather than breaking test collection at import time).
"""
import os

import pytest

_INDEX = "tests/benchmark/test_idx_only/genes_idx/ref.mk"
_OPTS = {"min_seed_len": 14, "mark_secondary": True, "softclip_supplementary": True}


@pytest.fixture(scope="module")
def bwa():
    if not os.path.exists(_INDEX):
        pytest.skip("fork BWA benchmark index not present")
    try:
        from bwamem import BwaAligner
    except ImportError:
        pytest.skip("bwamem not installed")
    return BwaAligner(_INDEX, **_OPTS)


def test_fork_aligns_some_reads(bwa):
    """The loaded aligner can align (reads) without error."""
    res = bwa.align_raw("ACGTACGTACGTACGT" * 5)
    assert isinstance(res, list)
