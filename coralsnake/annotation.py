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


def region_aligned_breaks(
    gene_splits: tuple[float, float, float], bin_number: int = 100
) -> np.ndarray:
    """
    Build bin edges aligned to the 5'UTR / CDS / 3'UTR region boundaries so that
    no single bin straddles a region junction.

    Feature positions are normalised per-region, so the region boundaries live at
    ``gene_splits[0]`` and ``gene_splits[0] + gene_splits[1]``. With the *uniform*
    ``linspace(0, 1, N+1)`` breaks one bin can span a junction, which pulls e.g.
    the first CDS-start sites into the last 5'UTR bin and inflates the boundary
    bin (the spurious start-codon "peak"). Placing the edges exactly on the
    boundaries keeps every bin wholly inside a single region.

    Within each region the bins are allocated proportionally to the region's
    normalized fraction, so every bin covers the same real fraction of the
    aligned transcript.

    Falls back to uniform breaks when the splits are all zero.
    """
    s0, s1, s2 = (float(x) for x in gene_splits)
    if s0 == 0 and s1 == 0 and s2 == 0:
        return np.linspace(0, 1, bin_number + 1)

    s0 = min(s0, 1.0)
    b2 = min(1.0, s0 + s1)
    if b2 < s0:
        b2 = s0

    if b2 >= 1.0:
        bounds = [0.0, 1.0] if s0 <= 0 else sorted({0.0, s0, 1.0})
    else:
        bounds = sorted({0.0, s0, b2, 1.0})

    if len(bounds) < 2:
        return np.linspace(0, 1, bin_number + 1)

    chunks: list[np.ndarray] = []
    for i in range(len(bounds) - 1):
        lo, hi = bounds[i], bounds[i + 1]
        width = hi - lo
        if width <= 0:
            continue
        n = max(1, int(round(width * bin_number)))
        seg = np.linspace(lo, hi, n + 1)
        if chunks:
            chunks.append(seg[1:])
        else:
            chunks.append(seg)

    breaks = np.concatenate(chunks)
    # Defensive: keep strictly increasing (guards against any rounding dup).
    return breaks[np.concatenate([np.array([True]), np.diff(breaks) > 0])]


def normalize_positions(
    annotated_sites: pl.DataFrame,
    split_strategy: str = "median",
    bin_number: int = 100,
    weight_col_index: list[int] | None = None,
    gene_splits: tuple | None = None,
    metric: str = "sum",
) -> tuple[pl.DataFrame, dict, tuple]:
    """
    Normalize transcript positions to relative feature positions (0-1 scale).
    Returns the normalized DataFrame and the gene splits.

    Bins are aligned to the 5'UTR / CDS / 3'UTR boundaries (see
    :func:`region_aligned_breaks`) so no bin straddles a region junction.

    ``metric`` selects the per-bin weight aggregation: ``"sum"`` yields
    ``count_<col>`` (total weighted signal, may be inflated in compressed
    short regions) and ``"mean"`` yields ``mean_<col>`` (signal per site,
    comparable across regions). Both are always emitted; ``metric`` only
    affects the returned value's naming preference for callers such as the
    CLI that export a single column.
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
                breaks=region_aligned_breaks(gene_splits, bin_number).tolist()
            )
        )
    )
    breaks = region_aligned_breaks(gene_splits, bin_number)
    n2c = {}
    if weight_col_index is None or len(weight_col_index) == 0:
        bin_counts, _ = np.histogram(
            gene_bins["feature_pos"],
            weights=gene_bins["feature_weight"],
            bins=breaks,
        )
        n2c["count"] = bin_counts
    else:
        # Raw site count per bin, used to derive the per-site mean (a region-
        # comparable metric that is not inflated by short-region compression).
        site_counts, _ = np.histogram(gene_bins["feature_pos"], bins=breaks)
        for col_index in weight_col_index:
            col_name = annotated_sites.columns[col_index]
            # The weight column may still be Unicode (a score column read as
            # text from the input); cast before multiplying, and treat
            # unparseable values as 0 so they do not crash the histogram.
            weight = gene_bins[col_name].cast(pl.Float64, strict=False).fill_null(0.0)
            bin_counts, _ = np.histogram(
                gene_bins["feature_pos"],
                weights=gene_bins["feature_weight"] * weight,
                bins=breaks,
            )
            n2c[f"count_{col_name}"] = bin_counts
            n2c[f"mean_{col_name}"] = np.divide(
                bin_counts,
                site_counts,
                out=np.zeros_like(bin_counts, dtype=float),
                where=site_counts > 0,
            )
    bin_midpoints = (breaks[:-1] + breaks[1:]) / 2
    gene_bins = pl.DataFrame({"feature_midpoint": bin_midpoints, **n2c})
    return gene_bins, gene_stats, gene_splits


def reference_point_positions(reference_df: pl.DataFrame, point: str) -> pl.DataFrame:
    """
    Build the per-transcript reference-feature positions for the reference-point
    metagene, in transcript-relative 5'->3' 0-based coordinates.

    Args:
        reference_df: the exon-level reference frame (from ``load_gtf`` /
            ``load_reference``), carrying per-exon ``transcript_id``,
            ``Start_exon``, ``transcript_length``, ``start_codon_pos`` and
            ``stop_codon_pos``.
        point: one of ``tss``, ``tes``, ``start_codon``, ``stop_codon``,
            ``last_exon_start``, ``exon_junction``.

    Returns:
        * single-point types: one row per transcript with an Int64 ``ref_pos``.
        * ``exon_junction``: one row per transcript with ``junctions`` = the
          sorted list of *internal* splice-junction positions (every exon 5'
          boundary except the transcript 5' end / TSS; the last exon's 5' boundary
          is the final splice junction). Single-exon transcripts have no internal
          junction and are omitted.
    """
    if point == "tss":
        return reference_df.group_by("transcript_id").agg(
            ref_pos=pl.lit(0, dtype=pl.Int64)
        )
    if point == "tes":
        return reference_df.group_by("transcript_id").agg(
            ref_pos=pl.col("transcript_length").first().cast(pl.Int64)
        )
    if point == "start_codon":
        return reference_df.group_by("transcript_id").agg(
            ref_pos=pl.col("start_codon_pos").first().cast(pl.Int64)
        )
    if point == "stop_codon":
        return reference_df.group_by("transcript_id").agg(
            ref_pos=pl.col("stop_codon_pos").first().cast(pl.Int64)
        )
    if point == "last_exon_start":
        # 5' boundary of the 3'-most exon (= last in 5'->3' order).
        return reference_df.group_by("transcript_id").agg(
            ref_pos=pl.col("Start_exon").max().cast(pl.Int64)
        )
    if point == "exon_junction":
        return (
            reference_df.group_by("transcript_id")
            .agg(_starts=pl.col("Start_exon").sort())
            .with_columns(junctions=pl.col("_starts").list.slice(1))
            .select("transcript_id", "junctions")
        )
    raise ValueError(f"Unknown reference point: {point}")


def _assign_nearest_junction(
    sites: pl.DataFrame, junctions_by_tx: pl.DataFrame
) -> pl.DataFrame:
    """Assign each mapped site to its nearest internal splice junction.

    Adds an Int64 ``ref_pos`` column = the transcript coordinate of the closest
    internal junction (by abs distance), vectorized per-transcript via
    ``searchsorted`` so we avoid a site x junction Cartesian blow-up. Sites whose
    transcript has no internal junction are dropped.

    Args:
        sites: annotated sites with ``transcript_id`` and ``transcript_pos``.
        junctions_by_tx: output of :func:`reference_point_positions` for
            ``point="exon_junction"`` (``transcript_id`` + ``junctions`` list).
    """
    jmap: dict = {}
    for row in junctions_by_tx.iter_rows(named=True):
        arr = np.sort(np.asarray(row["junctions"], dtype=np.int64))
        if arr.size:
            jmap[row["transcript_id"]] = arr

    parts = sites.partition_by("transcript_id")
    pieces = []
    for df in parts:
        tid = df["transcript_id"][0]
        jarr = jmap.get(tid)
        if jarr is None:
            continue
        pos = df["transcript_pos"].to_numpy().astype(np.int64)
        idx = np.clip(np.searchsorted(jarr, pos, side="left"), 0, len(jarr))
        lo = np.where(idx > 0, jarr[np.maximum(idx - 1, 0)], -(10**15))
        hi = np.where(idx < len(jarr), jarr[np.minimum(idx, len(jarr) - 1)], 10**15)
        best = np.where(np.abs(pos - lo) <= np.abs(pos - hi), lo, hi)
        pieces.append(df.with_columns(pl.Series("ref_pos", best, dtype=pl.Int64)))
    if not pieces:
        return sites.with_columns(pl.lit(None, dtype=pl.Int64).alias("ref_pos"))
    return pl.concat(pieces)


def normalize_point_positions(
    annotated_sites: pl.DataFrame,
    reference_points: pl.DataFrame,
    point: str,
    span_before: int,
    span_after: int,
    bin_number: int = 100,
    weight_col_index: list[int] | None = None,
    metric: str = "sum",
) -> tuple[pl.DataFrame, int]:
    """
    Reference-point metagene binning (DeepTools ``computeMatrix`` reference-point).

    Each mapped site is placed by its bp distance to the per-transcript reference
    feature (``transcript_pos - ref_pos``) and aggregated over
    ``[-span_before, +span_after]`` into ``bin_number`` uniform bins (region
    boundaries are not used here, so breaks are uniform over ``[0, 1]``). Sites
    whose distance is outside the span, or that have no reference feature, are
    dropped.

    Returns ``(gene_bins, n_dropped)``. ``gene_bins`` carries ``feature_midpoint``
    (normalized ``[0, 1]``), ``distance_bp`` (signed bp, 0 = reference feature)
    and ``count_*`` / ``mean_*`` columns, matching ``normalize_positions``.
    """
    sites = annotated_sites.with_columns(
        transcript_pos=(pl.col("transcript_start") + pl.col("transcript_end")) // 2
    ).filter(pl.col("transcript_pos").is_not_null())

    if point == "exon_junction":
        sites = _assign_nearest_junction(sites, reference_points)
    else:
        sites = sites.join(reference_points, on="transcript_id", how="inner")
        sites = sites.filter(pl.col("ref_pos").is_not_null())

    total_span = span_before + span_after
    n_total = sites.height
    sites = (
        sites.with_columns(distance=pl.col("transcript_pos") - pl.col("ref_pos"))
        .filter(
            (pl.col("distance") >= -span_before) & (pl.col("distance") <= span_after)
        )
        .with_columns(
            feature_pos=(pl.col("distance") + span_before) / total_span,
            feature_weight=pl.lit(1.0),
        )
    )
    n_dropped = n_total - sites.height

    breaks = np.linspace(0, 1, bin_number + 1)
    n2c: dict = {}
    if weight_col_index is None or len(weight_col_index) == 0:
        counts, _ = np.histogram(
            sites["feature_pos"], weights=sites["feature_weight"], bins=breaks
        )
        n2c["count"] = counts
    else:
        site_counts, _ = np.histogram(sites["feature_pos"], bins=breaks)
        for col_index in weight_col_index:
            col_name = annotated_sites.columns[col_index]
            weight = sites[col_name].cast(pl.Float64, strict=False).fill_null(0.0)
            counts, _ = np.histogram(
                sites["feature_pos"],
                weights=sites["feature_weight"] * weight,
                bins=breaks,
            )
            n2c[f"count_{col_name}"] = counts
            n2c[f"mean_{col_name}"] = np.divide(
                counts,
                site_counts,
                out=np.zeros_like(counts, dtype=float),
                where=site_counts > 0,
            )
    midpoints = (breaks[:-1] + breaks[1:]) / 2
    gene_bins = pl.DataFrame(
        {
            "feature_midpoint": midpoints,
            "distance_bp": (midpoints * total_span) - span_before,
            **n2c,
        }
    )
    return gene_bins, n_dropped


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
