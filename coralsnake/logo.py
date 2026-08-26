#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# DNA/RNA motif-logo plotting for the coralsnake CLI (`coralsnake logo`).
#
# Migrated from the standalone `motiflogo` package. matplotlib is an OPTIONAL
# dependency (extra: "coralsnake[plot]") and is imported lazily so that simply
# installing coralsnake does not pull in the heavy matplotlib stack.

from collections import defaultdict

import numpy as np

__all__ = ["Mlogo", "plot_motif_score"]

COLOR_SCHEME = {
    "A": "r",
    "C": "b",
    "G": "orange",
    "N": "black",
    "T": "darkgreen",
    "U": "darkgreen",
}


def _require_plotting():
    """Import matplotlib lazily, raising a helpful error if it is missing.

    matplotlib is intentionally an optional dependency (heavy). Users should
    ``pip install coralsnake[plot]`` to enable visualization.
    """
    try:
        import matplotlib as mpl
        import matplotlib.transforms as mtransforms
        from matplotlib.font_manager import FontProperties
        from matplotlib.patches import PathPatch
        from matplotlib.pyplot import gca
        from matplotlib.text import TextPath

        return mpl, mtransforms, FontProperties, PathPatch, gca, TextPath
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "Motif-logo plotting requires the optional 'matplotlib' dependency.\n"
            "Install it with:  pip install 'coralsnake[plot]'"
        ) from e


_PLOT_CTX = {}


def _plot_context():
    """Build (and cache) the matplotlib plotting helpers used by the logo.

    matplotlib is optional, so these are only constructed on first use.
    """
    if _PLOT_CTX:
        return _PLOT_CTX
    _, mtransforms, FontProperties, PathPatch, gca, TextPath = _require_plotting()
    fp = FontProperties(family="Arial", weight="bold")
    globscale = 1.375
    letters = {
        "A": TextPath((-0.350, 0), "A", size=1.015, prop=fp),
        "C": TextPath((-0.366, 0), "C", size=1.00, prop=fp),
        "G": TextPath((-0.384, 0), "G", size=1.00, prop=fp),
        "N": TextPath((-0.300, 0), "N", size=1.00, prop=fp),
        "T": TextPath((-0.305, 0), "T", size=1.00, prop=fp),
        "U": TextPath((-0.305, 0), "U", size=1.01, prop=fp),
    }
    _PLOT_CTX.update(
        {
            "mtransforms": mtransforms,
            "PathPatch": PathPatch,
            "gca": gca,
            "TextPath": TextPath,
            "FontProperties": FontProperties,
            "letters": letters,
            "globscale": globscale,
        }
    )
    return _PLOT_CTX


def letterAt(letter, x, y, yscale, color, ax, **kwargs):
    ctx = _plot_context()
    letters = ctx["letters"]
    text = letters.get(
        letter,
        ctx["TextPath"](
            (-0.305, 0),
            letter,
            size=1,
            prop=ctx["FontProperties"](family="Arial", weight="bold"),
        ),
    )
    t = (
        ctx["mtransforms"]
        .Affine2D()
        .scale(1 * ctx["globscale"], yscale * ctx["globscale"])
        + ctx["mtransforms"].Affine2D().translate(x, y)
        + ax.transData
    )
    p = ctx["PathPatch"](text, lw=0, fc=color, transform=t, **kwargs)
    if ax is not None:
        ax.add_artist(p)
    return p


def plot_motif_score(
    motif_score,
    ax=None,
    colors: dict = {},
    mask_index: int | list = [],
    mask_base: str | list = [],
    **kwargs,
):
    """Plot a computed motif score matrix as a sequence logo.

    Args:
        motif_score: list of lists of (base, score) tuples, one list per column.
        ax: optional matplotlib axes. If None, uses the current axes.
        colors: optional dict overriding the default letter colors.
        mask_index: single int or list of ints (column indices) to mask in grey.
        mask_base: single base str (mask that base everywhere) or list of bases.
    """
    _, _, _, _, gca, _ = _require_plotting()
    mask_color = "silver"
    if ax is None:
        ax = gca()
    # turn colors into a defaultdict, defaulting to the mask color
    c0 = defaultdict(lambda: mask_color) | COLOR_SCHEME
    if not isinstance(mask_index, list):
        mask_index = [mask_index]
    colors_dict = {
        idx: (c0 | colors) if idx not in mask_index else c0
        for idx in range(len(motif_score))
    }
    if isinstance(mask_base, str):
        for idx in colors_dict:
            colors_dict[idx][mask_base] = mask_color
    else:
        for idx, base in zip(colors_dict, list(mask_base)):
            colors_dict[idx][base] = mask_color

    x = 1
    maxi = 0
    for idx, scores in enumerate(motif_score):
        y = 0
        for base, score in scores:
            c = colors_dict[idx][base]
            letterAt(base, x, y, score, c, ax, **kwargs)
            y += score
        x += 1
        maxi = max(maxi, y)
    ax.set_xlim((0, x))
    ax.set_ylim((0, maxi))
    return ax


class Mlogo:
    """Build a sequence-logo score matrix from a collection of motifs.

    Mirrors the standalone `motiflogo.Mlogo` API.
    """

    def __init__(
        self,
        motifs: list = [],
        weights: list = [],
        scores: list = [],
        t2u=True,
        to2bit=True,
        normed=False,
    ):
        self.motifs = motifs

        # weight for each motif, same length as motif_list
        if len(weights) == 0:
            _weights = np.ones(len(motifs))
        else:
            _weights = np.nan_to_num(weights, nan=0)
            _weights = weights - np.min(_weights)
        weights = list(_weights)
        self.weights = weights

        self.t2u = t2u
        self.to2bit = to2bit
        self.normed = normed
        if len(motifs) > 0:
            self.scores = self._motif_to_score()
        else:
            # score is a list of list of tuple, check if the input is valid
            if all(isinstance(i, list) for i in scores):
                self.scores = scores
            else:
                raise ValueError("The input score is not valid")

    def _motif_to_score(self):
        motif_list = list(self.motifs)
        weights = self.weights
        t2u = self.t2u
        to2bit = self.to2bit
        normed = self.normed
        motif_len_all = list(set(len(m) for m in motif_list))
        motif_num = len(motif_list)
        if len(motif_len_all) == 1:
            motif_len = motif_len_all[0]
        else:
            raise ValueError("The motifs are not the same in length")

        if t2u:
            motif_list = [m.replace("T", "U") for m in motif_list]

        mm = []
        for i in range(motif_len):
            m = defaultdict(float)
            for j in range(motif_num):
                m[motif_list[j][i]] += weights[j]
            mm.append(sorted(m.items(), key=lambda x: x[1], reverse=True))

        if to2bit:
            mmm = []
            for m in mm:
                ss = (
                    2 + sum(i / motif_num * np.log2(i / motif_num) for _, i in m)
                ) / motif_num
                mmm.append([(b, i * ss) for b, i in m])
            return mmm

        if normed:
            total = sum(weights)
            mmm = []
            for m in mm:
                mmm.append([(b, i / total) for b, i in m])
            return mmm

        return mm

    def plot(
        self,
        ax=None,
        colors: dict = {},
        mask_index: int | list = [],
        mask_base: str | list = [],
        **kwargs,
    ):
        plot_motif_score(
            self.scores,
            ax=ax,
            colors=colors,
            mask_index=mask_index,
            mask_base=mask_base,
            **kwargs,
        )
        return ax
