#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regression tests for coralsnake.overlap (annotate_with_features / bins)."""

import polars as pl

from coralsnake.overlap import annotate_with_features, calculate_bin_statistics


def _data():
    inp = pl.DataFrame(
        {
            "Chromosome": ["chr1", "chr1"],
            "Start": [10, 20],
            "End": [30, 40],
            "Strand": ["+", "-"],
        }
    )
    feat = pl.DataFrame(
        {
            "Chromosome": ["chr1", "chr1"],
            "Start": [0, 0],
            "End": [50, 50],
            "Strand": ["+", "-"],
            "Transcript": ["T1", "T2"],
            "Type": ["CDS", "CDS"],
            "Name_b": ["g1", "g2"],
            "len_of_window": [50, 50],
            "len_of_feature": [50, 50],
            "frac_of_feature": [0.0, 0.0],
        }
    )
    return inp, feat


def test_annotate_with_features_no_strand():
    inp, feat = _data()
    df, score = annotate_with_features(
        inp, feat, bin_number=10, type_ratios=[1, 1, 1], by_strand=False
    )
    assert {"Chromosome", "Start", "End", "d_norm", "Strand"} <= set(df.columns)
    assert df.height == 2
    assert "count" in score.columns


def test_annotate_with_features_by_strand():
    # Regression: the by_strand branch previously referenced undefined variables.
    inp, feat = _data()
    df, _score = annotate_with_features(
        inp, feat, bin_number=10, type_ratios=[1, 1, 1], by_strand=True
    )
    assert df.height == 2


def test_annotate_with_features_annot_name():
    # annot_name=True should include the feature "Name" column in the output.
    inp, feat = _data()
    df, _ = annotate_with_features(
        inp, feat, bin_number=10, type_ratios=[1, 1, 1], annot_name=True
    )
    assert "Name" in df.columns


def test_calculate_bin_statistics():
    df = calculate_bin_statistics(
        [0.1, 0.5, 0.9], weights=[1.0, 2.0, 3.0], num_bins=10
    )
    assert {"count", "sum", "mean", "index"} <= set(df.columns)
    assert df["count"].sum() == 3.0
