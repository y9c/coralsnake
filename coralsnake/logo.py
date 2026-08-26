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


from .utils import require_plotting


def _require_plotting():
    """Import the matplotlib submodules needed to draw a sequence logo."""
    mpl = require_plotting(backend=None)
    import matplotlib.transforms as mtransforms
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch
    from matplotlib.pyplot import gca
    from matplotlib.text import TextPath

    return mpl, mtransforms, FontProperties, PathPatch, gca, TextPath


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
            # Min-shift from the CLEANED array, otherwise any NaN in the original
            # `weights` would propagate back in (NaN - min == NaN) and defeat
            # the nan_to_num above.
            _weights = _weights - np.min(_weights)
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
        motifs_list = list(self.motifs)
        weights = self.weights
        t2u = self.t2u
        to2bit = self.to2bit
        normed = self.normed

        motif_lens = {len(m) for m in motifs_list}
        if len(motif_lens) != 1:
            raise ValueError("The motifs are not the same in length")
        motif_len = motif_lens.pop()
        motif_num = len(motifs_list)

        if motif_len == 0 or motif_num == 0:
            return []

        if t2u:
            motifs_list = [m.replace("T", "U") for m in motifs_list]

        symbols, codes, weights_arr, first_index = self._encode(motifs_list, weights)
        nsym = len(symbols)

        # Vectorized per-column weighted counts via a single np.bincount
        # (far faster than np.add.at for large N). key = position*nsym + symbol.
        key = (codes + np.arange(motif_len)[None, :] * nsym).ravel()
        flat_w = np.broadcast_to(weights_arr[:, None], (motif_num, motif_len)).ravel()
        counts = np.bincount(key, weights=flat_w, minlength=motif_len * nsym).reshape(
            motif_len, nsym
        )

        # Per-column present symbols = those appearing in >=1 motif (weight may be 0).
        present_by_col = first_index < motif_num

        # Build per-column (base, score) lists with stable ordering:
        # primary = count descending; tie = first motif index that carried the base.
        mm = []
        for i in range(motif_len):
            present = np.flatnonzero(present_by_col[i])
            order = np.lexsort((first_index[i][present], -counts[i][present]))
            mm.append([(symbols[k], float(counts[i][k])) for k in present[order]])

        if to2bit:
            # Information (2-bit) correction per column, vectorized.
            p = counts / motif_num
            with np.errstate(divide="ignore", invalid="ignore"):
                entropy = np.where(p > 0, p * np.log2(p), 0.0)
            info = (2.0 + entropy.sum(axis=1)) / motif_num
            return [[(b, float(v) * s) for b, v in col] for col, s in zip(mm, info)]

        if normed:
            total = float(np.sum(weights_arr))
            if total == 0.0:
                # Original motiflogo would ZeroDivisionError here; return zeros
                # gracefully instead (all weights equal after min-shift).
                return [[(b, 0.0) for b, _ in col] for col in mm]
            return [[(b, v / total) for b, v in col] for col in mm]

        return mm

    @staticmethod
    def _encode(motifs: list[str], weights: list[float]):
        """Vectorize a list of motifs into (symbols, codes, weights, first_index).

        Returns
        -------
        symbols : list[str]
            Distinct bases, in sorted order (reproducible).
        codes : np.ndarray[int64], shape (N, L)
            Symbol code per motif/position (0..len(symbols)-1).
        weights : np.ndarray[float64], shape (N,)
            Per-motif weights.
        first_index : np.ndarray[int64], shape (L, len(symbols))
            First motif index that carries each symbol at each position, or
            ``N`` (sentinel) when a symbol is absent at that column.
        """
        n = len(motifs)
        motif_len = len(motifs[0])

        # Concatenate all motif strings and reinterpret as a byte buffer so the
        # (N, L) code matrix is produced by one vectorized gather (no list-of-
        # lists and no 2D np.unique).
        buf = "".join(motifs).encode("ascii", errors="ignore")
        raw = np.frombuffer(buf, dtype=np.uint8)  # shape (n * L,)
        uniq = np.unique(raw)
        symbols = [chr(int(x)) for x in uniq]
        table = np.zeros(256, dtype=np.int64)
        table[uniq] = np.arange(len(uniq))
        codes = table[raw].reshape(n, motif_len)

        weights_arr = np.asarray(weights, dtype=np.float64)
        if weights_arr.ndim == 0:
            weights_arr = np.full(n, float(weights_arr))

        # First motif index that carries each symbol at each column.
        first_index = np.full((motif_len, len(symbols)), n, dtype=np.int64)
        for c in range(len(symbols)):
            rows, cols = np.nonzero(codes == c)
            if rows.size:
                np.minimum.at(first_index[:, c], cols, rows)

        return symbols, codes, weights_arr, first_index

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
