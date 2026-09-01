#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Map to Local Coordinates - Convert global genomic coordinates to local reference coordinates

This module implements map_to_local functionality using ruranges,
inspired by pyranges1's map_to_local but using Rust-based ruranges operations.
"""

import numpy as np
import polars as pl
from ruranges.numpy import group_cumsum, overlaps

from .utils import interval_groups


def _strand_aware_cumsum(reference: pl.DataFrame, ref_id_col: str) -> pl.DataFrame:
    """Attach strand-aware cumulative exon offsets to a reference.

    Uses ruranges' vectorized ``group_cumsum`` (Rust) instead of a slow
    per-group Python apply. Returns the reference with ``_cumsum_start`` and
    ``_cumsum_end`` columns giving the 5'->3' running transcript offsets.

    Note: ruranges' ``negative_strand`` encodes *plus* strand as True.
    """
    starts = reference["Start"].cast(pl.Int64).to_numpy()
    ends = reference["End"].cast(pl.Int64).to_numpy()
    strands = reference["Strand"].to_numpy()

    # Build a stable numeric group id per ref_id_col value.
    uniq = reference[ref_id_col].unique().to_list()
    group_of = {v: i for i, v in enumerate(uniq)}
    groups = np.fromiter(
        (group_of[v] for v in reference[ref_id_col].to_list()),
        dtype=np.uint32,
        count=len(reference),
    )
    # ruranges convention: minus strand -> False, plus strand -> True.
    negative_strand = np.array(strands == "+", dtype=bool)

    idx, cumsum_start, cumsum_end = group_cumsum(
        starts=starts,
        ends=ends,
        negative_strand=negative_strand,
        groups=groups,
        sort=False,
    )

    # Scatter the traversal-ordered cumsum back onto original row positions.
    full_start = np.zeros(len(starts), dtype=np.int64)
    full_end = np.zeros(len(starts), dtype=np.int64)
    full_start[idx] = cumsum_start
    full_end[idx] = cumsum_end

    return reference.with_columns(
        pl.Series("_cumsum_start", full_start, dtype=pl.Int64),
        pl.Series("_cumsum_end", full_end, dtype=pl.Int64),
    )


def map_to_local(
    query: pl.DataFrame,
    reference: pl.DataFrame,
    ref_id_col: str = "transcript_id",
    match_by: list[str] | str | None = None,
    keep_global_chrom: bool = False,
    keep_global_loc: bool = False,
) -> pl.DataFrame:
    """
    Map genomic intervals from global coordinates to local reference coordinates.

    This function transforms query intervals to local coordinates within reference intervals,
    similar to pyranges1's map_to_local. It uses ruranges for efficient overlap detection.

    Args:
        query: Polars DataFrame with genomic intervals to map (must have Chromosome, Start, End, Strand)
        reference: Polars DataFrame with reference intervals (must have Chromosome, Start, End, Strand, and ref_id_col)
        ref_id_col: Column name in reference that identifies unique references (e.g., "transcript_id", "gene_id")
        match_by: Column name(s) to match between query and reference (optional)
        keep_global_chrom: If True, keep original global chromosome as "Chromosome_global"
        keep_global_loc: If True, keep original global Start/End as "Start_global"/"End_global"

    Returns:
        Polars DataFrame with local coordinates

    Example:
        >>> # Map SNPs to transcript coordinates
        >>> snps = pl.DataFrame({
        ...     "Chromosome": ["chr1", "chr1"],
        ...     "Start": [1000, 2000],
        ...     "End": [1001, 2001],
        ...     "Strand": ["+", "+"],
        ... })
        >>> transcripts = pl.DataFrame({
        ...     "Chromosome": ["chr1"],
        ...     "Start": [900],
        ...     "End": [2100],
        ...     "Strand": ["+"],
        ...     "transcript_id": ["TX1"],
        ... })
        >>> result = map_to_local(snps, transcripts, ref_id_col="transcript_id")
    """
    # Validate required columns
    required_query_cols = ["Chromosome", "Start", "End", "Strand"]
    required_ref_cols = ["Chromosome", "Start", "End", "Strand", ref_id_col]

    for col in required_query_cols:
        if col not in query.columns:
            raise ValueError(f"Query DataFrame must have '{col}' column")

    for col in required_ref_cols:
        if col not in reference.columns:
            raise ValueError(f"Reference DataFrame must have '{col}' column")

    # Store original query columns for reordering
    original_query_cols = query.columns

    # Handle match_by parameter
    if match_by is not None:
        if isinstance(match_by, str):
            match_by = [match_by]
    else:
        match_by = []

    # Calculate cumulative positions using ruranges' vectorized, strand-aware
    # group_cumsum (replaces the hand-rolled group_by().map_groups apply).
    # negative_strand encodes 'True -> plus strand' per ruranges convention.
    ref_cum = _strand_aware_cumsum(reference, ref_id_col)
    ref_indexed = ref_cum.with_row_index("_ref_idx")

    query_indexed = query.with_row_index("_query_idx")

    # Prepare arrays for overlap detection
    query_starts = query_indexed["Start"].cast(pl.Int64).to_numpy()
    query_ends = query_indexed["End"].cast(pl.Int64).to_numpy()
    query_chroms = query_indexed["Chromosome"].to_numpy()
    # query_strands not needed for overlap grouping

    ref_starts = ref_indexed["Start"].cast(pl.Int64).to_numpy()
    ref_ends = ref_indexed["End"].cast(pl.Int64).to_numpy()
    ref_chroms = ref_indexed["Chromosome"].to_numpy()
    # ref_strands not needed for overlap grouping

    # Create group IDs by chromosome (strand is handled after overlap).
    query_groups, ref_groups = interval_groups(query_chroms, ref_chroms)

    # Find overlaps
    idx_query, idx_ref = overlaps(
        starts=query_starts,
        ends=query_ends,
        starts2=ref_starts,
        ends2=ref_ends,
        groups=query_groups,
        groups2=ref_groups,
    )

    if len(idx_query) == 0:
        # No overlaps found - return an empty result with the full expected
        # schema (extra query columns included, matching the non-empty path)
        result_cols = ["Chromosome", "Start", "End", "Strand"]
        for col in original_query_cols:
            if col not in ["Chromosome", "Start", "End", "Strand"]:
                result_cols.append(col)
        if keep_global_chrom:
            result_cols.append("Chromosome_global")
        if keep_global_loc:
            result_cols.extend(["Start_global", "End_global", "Strand_global"])
        return pl.DataFrame(
            {
                col: (
                    pl.Series([], dtype=query[col].dtype)
                    if col in query.columns
                    else pl.Series([], dtype=pl.Utf8)
                )
                for col in result_cols
            }
        )

    # Build overlapping pairs dataframe
    overlaps_df = pl.DataFrame(
        {
            "_query_idx": idx_query,
            "_ref_idx": idx_ref,
        }
    )

    # Join with original dataframes
    result = overlaps_df.join(query_indexed, on="_query_idx").join(
        ref_indexed, on="_ref_idx", suffix="_ref"
    )

    # Filter by match_by columns if specified
    if match_by:
        for col in match_by:
            if col in result.columns and f"{col}_ref" in result.columns:
                result = result.filter(pl.col(col) == pl.col(f"{col}_ref"))

    # Calculate intersection of query and reference intervals
    result = result.with_columns(
        [
            pl.max_horizontal(pl.col("Start"), pl.col("Start_ref")).alias(
                "_intersect_start"
            ),
            pl.min_horizontal(pl.col("End"), pl.col("End_ref")).alias("_intersect_end"),
        ]
    )

    # Store global coordinates if requested
    if keep_global_chrom:
        result = result.with_columns(pl.col("Chromosome").alias("Chromosome_global"))

    if keep_global_loc:
        result = result.with_columns(
            [
                pl.col("Start").alias("Start_global"),
                pl.col("End").alias("End_global"),
                pl.col("Strand").alias("Strand_global"),
            ]
        )

    # Transform coordinates to local reference coordinates
    # Handle strand-aware transformation
    ref_is_minus = result["Strand_ref"] == "-"

    # For minus strand references, transform differently
    result = result.with_columns(
        [
            pl.when(ref_is_minus)
            .then(
                pl.col("End_ref") - pl.col("_intersect_end") + pl.col("_cumsum_start")
            )
            .otherwise(
                pl.col("_intersect_start")
                - pl.col("Start_ref")
                + pl.col("_cumsum_start")
            )
            .alias("_local_start"),
            pl.when(ref_is_minus)
            .then(
                pl.col("End_ref") - pl.col("_intersect_start") + pl.col("_cumsum_start")
            )
            .otherwise(
                pl.col("_intersect_end") - pl.col("Start_ref") + pl.col("_cumsum_start")
            )
            .alias("_local_end"),
        ]
    )

    # Transform strand: if query and ref have same strand, result is +, otherwise -
    result = result.with_columns(
        pl.when(pl.col("Strand") == pl.col("Strand_ref"))
        .then(pl.lit("+"))
        .otherwise(pl.lit("-"))
        .alias("_local_strand")
    )

    # Replace coordinates with local coordinates
    # ref_id_col should be in the reference dataframe and may or may not have _ref suffix
    # ref_id_col is present in joined frame as provided

    result = result.with_columns(
        [
            pl.col(ref_id_col).alias("Chromosome"),  # ref_id becomes new chromosome
            pl.col("_local_start").alias("Start"),
            pl.col("_local_end").alias("End"),
            pl.col("_local_strand").alias("Strand"),
        ]
    )

    # Select output columns
    output_cols = ["Chromosome", "Start", "End", "Strand"]

    # Add any additional query columns (excluding coordinate columns)
    for col in original_query_cols:
        if (
            col not in ["Chromosome", "Start", "End", "Strand"]
            and col in result.columns
        ):
            output_cols.append(col)

    if keep_global_chrom:
        output_cols.append("Chromosome_global")

    if keep_global_loc:
        output_cols.extend(["Start_global", "End_global", "Strand_global"])

    # Filter to available columns and return
    available_cols = [col for col in output_cols if col in result.columns]
    return result.select(available_cols)
