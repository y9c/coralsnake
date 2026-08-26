#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Regression tests for the performance-critical metagene internals.

These guard the vectorized rewrites (replacing slow group_by().map_groups()
Python applies with ruranges/polars vectorized ops) against silently returning
to a slow path or losing the strand-aware cumsum contract.
"""

import time

import numpy as np
import polars as pl
import pytest

pytest.importorskip(
    "coralsnake.metagene", reason="metagene deps not available", exc_type=ImportError
)

from coralsnake.metagene import load_gtf, map_to_local, map_to_transcripts
from coralsnake.metagene.map_to_local import _strand_aware_cumsum

DATA = "tests/data/R64-1-1.release57.gtf"


def _sites(tmp_path, n=5000, chrom="I"):
    """Build a sites DataFrame sampling exon midpoints on a chromosome."""
    ref = load_gtf(DATA)
    exons = ref.filter(pl.col("Chromosome") == chrom)
    mid = (exons["Start"].to_numpy() + exons["End"].to_numpy()) // 2
    idx = np.random.default_rng(0).integers(0, len(mid), size=n)
    return pl.DataFrame(
        {
            "Chromosome": [chrom] * n,
            "Start": mid[idx].tolist(),
            "End": (mid[idx] + 1).tolist(),
            "Strand": exons["Strand"].to_numpy()[idx].tolist(),
        }
    )


class TestVectorizedCumsum:
    def test_cumsum_is_strand_aware(self):
        """group_cumsum path must reproduce the hand-rolled strand-aware offsets."""
        ref = load_gtf(DATA).filter(pl.col("Chromosome") == "I")
        cum = _strand_aware_cumsum(ref, "transcript_id")
        assert "_cumsum_start" in cum.columns and "_cumsum_end" in cum.columns
        # Every transcript's first (5') exon starts at 0 cumsum.
        first = (
            cum.sort(["transcript_id", "_cumsum_start"])
            .group_by("transcript_id")
            .first()
        )
        assert (first["_cumsum_start"] == 0).all()

    def test_group_cumsum_exon_lengths_match(self):
        ref = load_gtf(DATA).filter(pl.col("Chromosome") == "I")
        cum = _strand_aware_cumsum(ref, "transcript_id")
        # cumulative end - cumulative start == exon length for every row.
        lengths = (cum["End"] - cum["Start"]).to_list()
        offsets = (cum["_cumsum_end"] - cum["_cumsum_start"]).to_list()
        assert lengths == offsets


class TestPerformance:
    def test_map_to_transcripts_stays_fast(self, tmp_path):
        """Guard against regressing to the slow group_by().map_groups apply.
        The vectorized rewrite is ~20x faster; keep a generous ceiling."""
        ref = load_gtf(DATA)
        sites = _sites(tmp_path, n=5000)
        t0 = time.perf_counter()
        annot = map_to_transcripts(sites, ref)
        elapsed = time.perf_counter() - t0
        assert annot.height == sites.height
        # Ceiling: old python-apply path took ~250ms for 5k sites;
        # vectorized path is ~10ms. Allow 20x headroom to avoid flakiness.
        assert elapsed < 2.0, f"map_to_transcripts too slow ({elapsed:.3f}s)"

    def test_map_to_local_stays_fast(self, tmp_path):
        ref = load_gtf(DATA)
        sites = _sites(tmp_path, n=5000)
        t0 = time.perf_counter()
        result = map_to_local(sites, ref, ref_id_col="transcript_id")
        elapsed = time.perf_counter() - t0
        assert result.height > 0
        assert elapsed < 1.0, f"map_to_local too slow ({elapsed:.3f}s)"
