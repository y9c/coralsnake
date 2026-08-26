#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test cases for the migrated map_to_local function.

Ported from the standalone `metagene` package into coralsnake
(coralsnake.map_to_local). Uses the high-performance
ruranges + polars stack.
"""

import polars as pl
import pytest

pytest.importorskip(
    "coralsnake.map_to_local", reason="metagene deps not available", exc_type=ImportError
)

from coralsnake.map_to_local import map_to_local  # noqa: E402


class TestMapToLocal:
    def test_basic_mapping_plus_strand(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1", "chr1", "chr1"],
                "Start": [100, 150, 250],
                "End": [110, 160, 260],
                "Strand": ["+", "+", "+"],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [90],
                "End": [300],
                "Strand": ["+"],
                "transcript_id": ["TX1"],
            }
        )
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        assert result.height == 3
        assert result["Chromosome"].to_list() == ["TX1", "TX1", "TX1"]
        assert result["Start"][0] == 10
        assert result["End"][0] == 20
        assert result["Start"][1] == 60
        assert result["End"][1] == 70

    def test_basic_mapping_minus_strand(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [100, 150],
                "End": [110, 160],
                "Strand": ["-", "-"],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [90],
                "End": [300],
                "Strand": ["-"],
                "transcript_id": ["TX1"],
            }
        )
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        assert result.height == 2
        assert result["Chromosome"].to_list() == ["TX1", "TX1"]
        assert result["Start"][0] == 190
        assert result["End"][0] == 200

    def test_multi_exon_transcript(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1", "chr1", "chr1"],
                "Start": [100, 250, 450],
                "End": [110, 260, 460],
                "Strand": ["+", "+", "+"],
                "query_id": ["Q1", "Q2", "Q3"],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [90, 240],
                "End": [200, 500],
                "Strand": ["+", "+"],
                "transcript_id": ["TX1", "TX1"],
            }
        )
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        assert result.height == 3
        assert all(result["Chromosome"] == "TX1")
        result = result.sort("query_id")
        assert result.filter(pl.col("query_id") == "Q1")["Start"][0] == 10
        assert result.filter(pl.col("query_id") == "Q2")["Start"][0] == 120
        assert result.filter(pl.col("query_id") == "Q3")["Start"][0] == 320

    def test_strand_transformation(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [100],
                "End": [110],
                "Strand": ["+"],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [90],
                "End": [300],
                "Strand": ["-"],
                "transcript_id": ["TX1"],
            }
        )
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        assert result["Strand"][0] == "-"

    def test_no_overlap(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [1000],
                "End": [1100],
                "Strand": ["+"],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [90],
                "End": [300],
                "Strand": ["+"],
                "transcript_id": ["TX1"],
            }
        )
        result = map_to_local(query, reference, ref_id_col="transcript_id")
        assert result.height == 0

    def test_keep_global_coordinates(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [100],
                "End": [110],
                "Strand": ["+"],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [90],
                "End": [300],
                "Strand": ["+"],
                "transcript_id": ["TX1"],
            }
        )
        result = map_to_local(
            query,
            reference,
            ref_id_col="transcript_id",
            keep_global_chrom=True,
            keep_global_loc=True,
        )
        assert "Chromosome_global" in result.columns
        assert result["Chromosome_global"][0] == "chr1"
        assert result["Start_global"][0] == 100

    def test_match_by_filter(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [100, 150],
                "End": [110, 160],
                "Strand": ["+", "+"],
                "gene_id": ["G1", "G2"],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [90, 140],
                "End": [200, 250],
                "Strand": ["+", "+"],
                "transcript_id": ["TX1", "TX2"],
                "gene_id": ["G1", "G2"],
            }
        )
        result = map_to_local(
            query, reference, ref_id_col="transcript_id", match_by="gene_id"
        )
        assert result.height == 2

    def test_missing_required_columns_error(self):
        query = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [100],
                "End": [110],
            }
        )
        reference = pl.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [90],
                "End": [300],
                "Strand": ["+"],
                "transcript_id": ["TX1"],
            }
        )
        with pytest.raises(
            ValueError, match="Query DataFrame must have 'Strand' column"
        ):
            map_to_local(query, reference, ref_id_col="transcript_id")
