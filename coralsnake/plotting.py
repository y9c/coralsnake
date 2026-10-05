#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Plotting Module - Visualization functions for metagene analysis.
#
# NOTE: matplotlib is an OPTIONAL dependency (extra: "coralsnake[plot]").
# It is imported lazily so that the core metagene analysis keeps working
# without the heavy matplotlib stack installed.

import polars as pl

from .utils import require_plotting, setup_logger

# Set up logger
logger = setup_logger(__name__)


def _require_plotting():
    """Set the non-interactive backend and import pyplot lazily."""
    require_plotting(backend="Agg")
    import matplotlib.pyplot as plt

    return plt


def plot_profile(
    gene_bins: pl.DataFrame,
    gene_splits: tuple[float, float, float],
    output_path: str,
    figsize: tuple[float, float] = (6.4, 4.2),
    gene_model: bool = True,
):
    """
    Create the metagene profile plot using matplotlib (optional dependency).

    A single axes holds the region-shaded metagene curve (x-axis at y=0). The
    gene-feature schematic (whole gene = narrow bar, CDS = wider bar, CDS box
    top attached to y=0, whole-gene bar centred on the CDS midline) is drawn in
    the figure margin directly below the x-axis, with the x tick labels and
    the axis label placed below it. Set ``gene_model=False`` to omit it.

    Args:
        gene_bins: DataFrame with ``feature_midpoint`` and ``count_*`` columns.
        gene_splits: Tuple of (5'UTR, CDS, 3'UTR) normalized region ratios.
        output_path: Destination file (PNG/SVG/PDF are supported).
        figsize: Figure size in inches.
        gene_model: When True, draw the gene-feature schematic below the x-axis.
    """
    plt = _require_plotting()
    from matplotlib.patches import Rectangle

    b1 = float(gene_splits[0])
    b2 = b1 + float(gene_splits[1])

    fig, ax = plt.subplots(figsize=figsize)  # one panel only

    # Region shading (5'UTR / CDS / 3'UTR)
    ax.axvspan(0, b1, color="#bbdefb", alpha=0.5, lw=0)
    ax.axvspan(b1, b2, color="#c8e6c9", alpha=0.5, lw=0)
    ax.axvspan(b2, 1, color="#ffe0b2", alpha=0.5, lw=0)

    # Profile curve(s)
    count_cols = [c for c in gene_bins.columns if c.startswith("count")]
    for col in count_cols:
        x = gene_bins["feature_midpoint"]
        y = gene_bins[col]
        ax.plot(x, y, linewidth=2.0, label=col.removeprefix("count_"))
        if len(count_cols) == 1:
            ax.fill_between(x, 0, y, alpha=0.22)
    if len(count_cols) == 1:
        ax.legend(frameon=False, loc="upper right")

    # Region boundaries
    for b in (b1, b2):
        ax.axvline(b, color="k", ls="--", lw=1.0, alpha=0.6)

    # Region labels on the curve
    ymax = ax.get_ylim()[1]
    lab = dict(
        ha="center",
        va="top",
        fontweight="bold",
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", pad=0.5),
    )
    ax.text(b1 / 2, ymax * 0.98, "5'UTR", **lab)
    ax.text((b2 + 1) / 2, ymax * 0.98, "3'UTR", **lab)

    ax.set_xlim(0, 1)
    ax.set_ylabel("Density", fontsize=11)
    ax.grid(alpha=0.2)

    # Layout in figure fractions; leave a band below the axes for the gene model
    # + tick labels + axis label. The x-axis (data y=0) sits at figure-y = BOT.
    L, R, TOP, BOT = 0.11, 0.97, 0.92, 0.30
    fig.subplots_adjust(left=L, right=R, top=TOP, bottom=BOT)
    W = R - L
    fx = lambda v: L + v * W  # data-x (0..1) -> figure-x

    if gene_model:
        # Hide the axis tick labels (we draw them, plus the gene model, below),
        # but keep the x-axis line at y=0.
        ax.tick_params(axis="x", labelbottom=False)
        y0 = BOT
        h = 0.05  # CDS box height (figure fractions), top attached to y0
        mid = y0 - h / 2
        # CDS box (top at y=0) and centred whole-gene bar, drawn in figure coords
        fig.add_artist(
            Rectangle(
                (fx(b1), y0 - h),
                fx(b2) - fx(b1),
                h,
                transform=fig.transFigure,
                color="black",
            )
        )
        fig.add_artist(
            Rectangle(
                (fx(0), mid - 0.008),
                fx(1) - fx(0),
                0.016,
                transform=fig.transFigure,
                color="black",
            )
        )
        fig.text(
            fx((b1 + b2) / 2),
            mid,
            "CDS",
            color="white",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
        )
        # x tick labels + axis label below the gene model
        for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
            fig.text(
                fx(t), y0 - h - 0.045, f"{t:.1f}", ha="center", va="top", fontsize=9
            )
        fig.text(
            fx(0.5), y0 - h - 0.11, "Normalized Gene Position", ha="center", fontsize=10
        )
    else:
        ax.set_xlabel("Normalized Gene Position", fontsize=10)

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
