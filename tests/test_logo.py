#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for the migrated motif-logo module (coralsnake.logo).

The score-matrix build (Mlogo) is pure numpy and always works. The actual
rendering requires the optional 'matplotlib' extra; those tests are skipped
when matplotlib is not installed.
"""

import numpy as np
import pytest

from coralsnake.logo import COLOR_SCHEME, Mlogo


class TestMlogo:
    def test_motif_to_scores_2bit(self):
        m = Mlogo(motifs=["ACGT", "ACGG"], to2bit=True, t2u=False)
        assert len(m.scores) == 4  # one score list per column
        # Each column has at least one (base, score) tuple
        assert all(len(col) >= 1 for col in m.scores)

    def test_t2u_conversion(self):
        # T should be converted to U in the generated logo
        m = Mlogo(motifs=["ATGC"], t2u=True, to2bit=False)
        bases = {b for col in m.scores for b, _ in col}
        assert "U" in bases
        assert "T" not in bases

    def test_normed_scores_sum_to_one(self):
        m = Mlogo(motifs=["ACGT", "ACGT", "CCGT"], normed=True, to2bit=False)
        for col in m.scores:
            total = sum(s for _, s in col)
            assert np.isclose(total, 1.0)

    def test_unequal_length_raises(self):
        import pytest

        with pytest.raises(ValueError, match="not the same in length"):
            Mlogo(motifs=["ACGT", "ACG"])

    def test_color_scheme(self):
        for base in ["A", "C", "G", "T", "U"]:
            assert base in COLOR_SCHEME


class TestLogoPlot:
    """Plotting requires the optional matplotlib extra."""

    def test_plot_requires_optional_dep(self):
        try:
            from coralsnake.logo import _require_plotting

            _require_plotting()
        except ImportError:
            pytest.skip("matplotlib not installed (optional 'plot' extra)")

    def test_plot_produces_file(self, tmp_path):
        from coralsnake.logo import _require_plotting

        try:
            _require_plotting()
        except ImportError:
            pytest.skip("matplotlib not installed (optional 'plot' extra)")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        out = tmp_path / "logo.png"
        m = Mlogo(motifs=["ACGT", "ACGG", "CCGT"])
        fig = plt.figure(figsize=(0.75 * len(m.scores), 2.5))
        ax = fig.gca()
        m.plot(ax=ax)
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        assert out.exists() and out.stat().st_size > 0
