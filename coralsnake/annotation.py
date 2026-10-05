#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Annotation and Metagene Analysis Functions

import numpy as np
import polars as pl
from ruranges.numpy import overlaps

from .utils import interval_groups


def map_to_transcripts(
    input_sites: pl.DataFrame, exon_ref: pl.DataFrame
) -> pl.DataFrame:
    """
    Annotate input sites with transcript information using exon reference.
    Returns a Polars DataFrame with transcript mapping.
    """
    # Add row index to input sites for later joining
    input_sites_indexed = input_sites.with_row_index("_tmp_row_index")

    # ruranges requires integer coordinate arrays. A nullable Int column (which
    # arises when a site coordinate failed to parse) is realised by numpy as
    # float64-with-NaN, which ruranges rejects. Drop sites without a finite
    # integer position from the overlap search, but keep their row index so the
    # final left-join still returns them (as unmapped/null annotations).
    # ``_pos`` is the positional index the overlap search returns (``input_idx``
    # indexes this filtered frame); ``_tmp_row_index`` is carried through to the
    # final left-join against the full input.
    sites_ok = input_sites_indexed.filter(
        pl.col("Start").is_not_null() & pl.col("End").is_not_null()
    ).with_row_index("_pos")

    # The best transcript per gene (min transcript_level, then max
    # transcript_length, then min transcript_id) is a property of the
    # reference: those three attributes are constant across a transcript's
    # exons. Resolve it ONCE, then restrict the overlap search itself to the
    # best transcripts' exons. That avoids the site x isoform blow-up in the
    # overlap output and the joins (the old code sorted the whole expanded
    # frame to pick a best transcript, and its peak memory was the bottleneck).
    # Verified equivalent: a site maps to the same (gene, transcript) set either
    # way.
    best_tx = (
        exon_ref.sort(
            ["gene_id", "transcript_level", "transcript_length", "transcript_id"],
            descending=[False, False, True, False],
        )
        .group_by("gene_id", maintain_order=True)
        .first()
        .select(["gene_id", "transcript_id"])
    )
    exon_ref_best = exon_ref.filter(
        pl.col("transcript_id").is_in(best_tx["transcript_id"].to_list())
    )
    exon_indexed = exon_ref_best.with_row_index("exon_idx")

    # Prepare arrays for overlap detection
    input_starts = sites_ok["Start"].cast(pl.Int64).to_numpy()
    input_ends = sites_ok["End"].cast(pl.Int64).to_numpy()
    input_chroms = sites_ok["Chromosome"].to_numpy()
    input_strands = sites_ok["Strand"].to_numpy()

    exon_starts = exon_indexed["Start"].cast(pl.Int64).to_numpy()
    exon_ends = exon_indexed["End"].cast(pl.Int64).to_numpy()
    exon_chroms = exon_indexed["Chromosome"].to_numpy()
    exon_strands = exon_indexed["Strand"].to_numpy()

    # Create group IDs combining chromosome and strand for strand-aware overlaps
    input_labels = np.char.add(input_chroms.astype(str), input_strands.astype(str))
    exon_labels = np.char.add(exon_chroms.astype(str), exon_strands.astype(str))
    input_groups, exon_groups = interval_groups(input_labels, exon_labels)

    # Find overlaps
    idx_exon, idx_input = overlaps(
        starts=exon_starts,
        ends=exon_ends,
        starts2=input_starts,
        ends2=input_ends,
        groups=exon_groups,
        groups2=input_groups,
    )

    # No overlaps is a valid outcome (e.g. every site is intergenic, or on a
    # contig absent from the reference): those sites come back as unmapped/null
    # annotations via the left-join below, matching the function's contract.
    # Only an empty reference itself is a configuration error.
    if len(exon_indexed) == 0:
        raise ValueError(
            "Reference contains no exons; cannot map any input sites. "
            "Check the reference/GTF."
        )

    # Build overlapping pairs dataframe
    overlaps_df = pl.DataFrame(
        {
            "exon_idx": idx_exon,
            "input_idx": idx_input,
        }
    )

    # Join with original dataframes. exon_indexed already contains only the
    # best transcripts' exons, so no further best-transcript filter is needed.
    annot = overlaps_df.join(exon_indexed, on="exon_idx", how="inner").join(
        sites_ok.select(
            ["_pos", "_tmp_row_index", "Chromosome", "Start", "End", "Strand"]
        ),
        left_on="input_idx",
        right_on="_pos",
        suffix="_qry",
    )

    # Add reference columns
    annot = annot.with_columns(
        pl.col("Chromosome").alias("Chromosome_ref"),
        pl.col("Start").alias("Start_ref"),
        pl.col("End").alias("End_ref"),
        pl.col("Strand").alias("Strand_ref"),
    )

    # Calculate new Start/End based on strand
    annot = annot.with_columns(
        [
            pl.when(pl.col("Strand_ref") == "+")
            .then(
                (pl.col("Start_qry") - pl.col("Start_ref")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .otherwise(
                (pl.col("End_ref") - pl.col("End_qry")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .alias("transcript_start"),
            pl.when(pl.col("Strand_ref") == "+")
            .then(
                (pl.col("End_qry") - pl.col("Start_ref")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .otherwise(
                (pl.col("End_ref") - pl.col("Start_qry")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .alias("transcript_end"),
        ]
    )

    annotation_cols = [
        "gene_id",
        "transcript_id",
        "transcript_start",
        "transcript_end",
        "transcript_length",
        "start_codon_pos",
        "stop_codon_pos",
        "exon_number",
        "Start_exon",
        "End_exon",
    ]
    annot = annot.select(["_tmp_row_index"] + annotation_cols)

    # Join annotation back to input_sites
    annotated_sites = (
        input_sites.with_row_index("_tmp_row_index")
        .drop(annotation_cols, strict=False)
        .join(annot, on="_tmp_row_index", how="left")
        .with_columns(
            pl.col("transcript_start").cast(pl.Int64),
            pl.col("transcript_end").cast(pl.Int64),
            pl.col("transcript_length").cast(pl.Int64),
            pl.col("start_codon_pos").cast(pl.Int64),
            pl.col("stop_codon_pos").cast(pl.Int64),
            pl.col("Start_exon").cast(pl.Int64),
            pl.col("End_exon").cast(pl.Int64),
        )
        .with_columns(record_id=pl.col("_tmp_row_index"))
        .drop("_tmp_row_index")
    )
    return annotated_sites


def calculate_gene_splits(
    annotated_sites: pl.DataFrame, split_strategy: str = "mean"
) -> tuple:
    """
    Calculate gene region splits (5'UTR, CDS, 3'UTR) from annotated sites.
    """
    df = (
        annotated_sites.select(
            "transcript_id", "transcript_length", "start_codon_pos", "stop_codon_pos"
        )
        .drop_nulls()
        .unique()
    )
    if split_strategy == "mean":
        # Cast to numeric first to ensure we get numeric types
        start_mean = df.select(pl.col("start_codon_pos").cast(pl.Float64).mean()).item()
        stop_mean = df.select(pl.col("stop_codon_pos").cast(pl.Float64).mean()).item()
        length_mean = df.select(
            pl.col("transcript_length").cast(pl.Float64).mean()
        ).item()

        len_5utr = float(start_mean or 0)
        len_cds = float(stop_mean or 0) + 3 - float(start_mean or 0)
        len_3utr = float(length_mean or 0) - (float(stop_mean or 0) + 3)
    elif split_strategy == "median":
        # Cast to numeric first to ensure we get numeric types
        start_median = df.select(
            pl.col("start_codon_pos").cast(pl.Float64).median()
        ).item()
        stop_median = df.select(
            pl.col("stop_codon_pos").cast(pl.Float64).median()
        ).item()
        length_median = df.select(
            pl.col("transcript_length").cast(pl.Float64).median()
        ).item()

        len_5utr = float(start_median or 0)
        len_cds = float(stop_median or 0) + 3 - float(start_median or 0)
        len_3utr = float(length_median or 0) - (float(stop_median or 0) + 3)
    else:
        raise ValueError(f"Unknown split_strategy: {split_strategy}")

    len_total = len_5utr + len_cds + len_3utr
    if len_total == 0:
        return (0.0, 0.0, 0.0)
    return len_5utr / len_total, len_cds / len_total, len_3utr / len_total


def normalize_positions(
    annotated_sites: pl.DataFrame,
    split_strategy: str = "median",
    bin_number: int = 100,
    weight_col_index: list[int] | None = None,
    gene_splits: tuple | None = None,
) -> tuple[pl.DataFrame, dict, tuple]:
    """
    Normalize transcript positions to relative feature positions (0-1 scale).
    Returns the normalized DataFrame and the gene splits.
    """
    # check if the "transcript_id", "transcript_start" and  "transcript_end" in the dataframe columns
    # use the mid point of transcript_start and transcript_end as transcript_pos

    annotated_sites_bins = (
        annotated_sites.with_columns(
            transcript_pos=(pl.col("transcript_start") + pl.col("transcript_end")) // 2
        )
        # record_id is a unique row index, so each site has weight 1.0;
        # the previous `1 / pl.len().over("record_id")` was a no-op window
        # pass over the whole frame.
        .with_columns(feature_weight=pl.lit(1.0))
        .with_columns(
            feature_type=pl.when(pl.col("transcript_pos").is_null())
            .then(pl.lit("None"))
            # noncoding transcript (no CDS start/stop): exclude from 5'UTR/CDS/3'UTR
            .when(
                pl.col("start_codon_pos").is_null() | pl.col("stop_codon_pos").is_null()
            )
            .then(pl.lit("None"))
            .when(pl.col("transcript_pos") < pl.col("start_codon_pos"))
            .then(pl.lit("5UTR"))
            # stop codon (stop_codon_pos .. stop_codon_pos+2) is part of the CDS;
            # 3'UTR starts at the first base after it (matches effect.py/_classify_exonic)
            .when(pl.col("transcript_pos") >= pl.col("stop_codon_pos") + 3)
            .then(pl.lit("3UTR"))
            .otherwise(pl.lit("CDS"))
        )
    )

    gene_stats = annotated_sites_bins.group_by("feature_type").agg(
        count=pl.col("feature_weight").sum()
    )
    gene_stats = dict(zip(gene_stats["feature_type"], gene_stats["count"]))

    # Allow callers to supply precomputed splits (e.g. computed on the full
    # gene population) so a region-filtered frame does not recompute the
    # 5'UTR/CDS/3'UTR boundaries from only the selected subset.
    if gene_splits is None:
        gene_splits = calculate_gene_splits(annotated_sites, split_strategy)

    gene_bins = (
        annotated_sites.with_columns(
            transcript_pos=(pl.col("transcript_start") + pl.col("transcript_end")) // 2
        )
        .filter(
            pl.col("transcript_pos").is_not_null()
            & pl.col("start_codon_pos").is_not_null()
            & pl.col("stop_codon_pos").is_not_null()
        )
        .with_columns(
            feature_pos=pl.when(pl.col("transcript_pos") < pl.col("start_codon_pos"))
            .then(
                pl.col("transcript_pos")
                / pl.col("start_codon_pos").clip(1, None)
                * gene_splits[0]
            )
            # CDS = [start_codon_pos, stop_codon_pos + 3) (stop codon is CDS)
            .when(pl.col("transcript_pos") >= pl.col("stop_codon_pos") + 3)
            .then(
                gene_splits[0]
                + gene_splits[1]
                + (pl.col("transcript_pos") - (pl.col("stop_codon_pos") + 3))
                / (pl.col("transcript_length") - (pl.col("stop_codon_pos") + 3)).clip(
                    1, None
                )
                * gene_splits[2]
            )
            .otherwise(
                gene_splits[0]
                + (pl.col("transcript_pos") - pl.col("start_codon_pos"))
                / (pl.col("stop_codon_pos") + 3 - pl.col("start_codon_pos")).clip(
                    1, None
                )
                * gene_splits[1]
            )
        )
        .with_columns(feature_weight=pl.lit(1.0))
        .with_columns(
            feature_bin=pl.col("feature_pos").cut(
                breaks=np.linspace(0, 1, bin_number + 1).tolist()
            )
        )
    )
    n2c = {}
    if weight_col_index is None or len(weight_col_index) == 0:
        bin_counts, _ = np.histogram(
            gene_bins["feature_pos"],
            weights=gene_bins["feature_weight"],
            bins=np.linspace(0, 1, bin_number + 1),
        )
        n2c["count"] = bin_counts
    else:
        for col_index in weight_col_index:
            col_name = annotated_sites.columns[col_index]
            # The weight column may still be Unicode (a score column read as
            # text from the input); cast before multiplying, and treat
            # unparseable values as 0 so they do not crash the histogram.
            weight = gene_bins[col_name].cast(pl.Float64, strict=False).fill_null(0.0)
            bin_counts, _ = np.histogram(
                gene_bins["feature_pos"],
                weights=gene_bins["feature_weight"] * weight,
                bins=np.linspace(0, 1, bin_number + 1),
            )
            n2c[f"count_{col_name}"] = bin_counts
    bin_midpoints = np.linspace(0, 1, bin_number + 1)[:-1] + 0.5 / bin_number
    gene_bins = pl.DataFrame({"feature_midpoint": bin_midpoints, **n2c})
    return gene_bins, gene_stats, gene_splits


def show_summary_stats(df: pl.DataFrame) -> str:
    """
    Generate summary statistics of the analysis.

    Args:
        df_normalized: Final DataFrame with all annotations

    Returns:
        A formatted string containing the summary statistics
    """
    # filter record with feature_type is not null, and show the proportion passed the filter
    total_passed = df.height
    total_sites = df.height
    pass_percentage = (total_passed / total_sites * 100) if total_sites > 0 else 0

    # Count by feature type
    feature_counts = df.group_by("feature_type").len().sort("feature_type")

    # Build feature distribution string
    feature_dist = []
    for row in feature_counts.iter_rows():
        feature_type, count = row
        percentage = (count / total_passed * 100) if total_passed > 0 else 0
        feature_dist.append(f"{feature_type}: {count} sites ({percentage:.1f}%)")

    # Calculate position statistics
    feature_positions = df["feature_pos"]
    pos_stats = []
    if len(feature_positions) > 0:
        pos_stats = [
            f"Mean: {feature_positions.mean():.3f}",
            f"Median: {feature_positions.median():.3f}",
            f"Min: {feature_positions.min():.3f}",
            f"Max: {feature_positions.max():.3f}",
        ]
    else:
        pos_stats = ["No valid position statistics available (all values are null)"]

    # Combine all parts into a single string
    summary = (
        f"Total sites passed the filter: {total_passed} / {total_sites} ({pass_percentage:.1f}%)\n\n"
        f"Feature Distribution:\n  "
        + "\n  ".join(feature_dist)
        + "\n\n"
        + "Position Statistics:\n  "
        + "\n  ".join(pos_stats)
    )

    return summary
