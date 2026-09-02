#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Tests for the machine-readable matrix exports added to the CLI:
#   * `coralsnake logo --matrix FILE`
#   * `coralsnake metagene --export-profile FILE`
#
# The logo matrix export must work WITHOUT matplotlib (that is the point of
# the flag), so the plotting extra is deliberately not required here.

import pytest

from coralsnake.logo import Mlogo

pytest.importorskip("coralsnake.cli", reason="cli deps not available", exc_type=ImportError)

from click.testing import CliRunner  # noqa: E402

from coralsnake.cli import cli  # noqa: E402

_GTF = """1\tt\tgene\t100\t400\t.\t+\t.\tgene_id "G1";
1\tt\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";
1\tt\texon\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";
1\tt\tstart_codon\t100\t102\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";
1\tt\tstop_codon\t398\t400\t.\t+\t.\tgene_id "G1"; transcript_id "G1.1";
"""


def _write_gtf(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(_GTF)
    return gtf


class TestLogoMatrix:
    def test_matrix_written_without_matplotlib(self, tmp_path):
        """--matrix works even when the plotting extra is missing."""
        out = tmp_path / "logo_matrix.tsv"
        res = CliRunner().invoke(
            cli, ["logo", "-m", "ACGT", "-m", "ACGU", "--matrix", str(out)]
        )
        assert res.exit_code == 0, res.output
        assert out.exists()
        header = out.read_text().splitlines()[0].split("\t")
        assert header[0] == "position"
        assert header[1:] == ["A", "C", "G", "T", "U"]

    def test_matrix_rows_match_motif_length(self, tmp_path):
        out = tmp_path / "logo_matrix.tsv"
        res = CliRunner().invoke(
            cli, ["logo", "-m", "ACGU", "--matrix", str(out)]
        )
        assert res.exit_code == 0, res.output
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 1 + 4  # header + one row per position
        assert [ln.split("\t")[0] for ln in lines[1:]] == [str(i) for i in range(1, 5)]

    def test_matrix_values_match_mlogo_scores(self, tmp_path):
        """The exported matrix equals Mlogo.scores, zero-padded per position."""
        motifs = ["ACGU", "CGUA"]
        # Build Mlogo the same way the CLI does (explicit weights -> min-shift)
        m = Mlogo(motifs=list(motifs), weights=[1.0, 1.0], t2u=True, to2bit=False)
        motifs_tsv = tmp_path / "motifs.tsv"
        motifs_tsv.write_text("".join(f"{mo}\t1\n" for mo in motifs))
        out = tmp_path / "logo_matrix.tsv"
        res = CliRunner().invoke(
            cli, ["logo", "-i", str(motifs_tsv), "--no-2bit", "--matrix", str(out)]
        )
        assert res.exit_code == 0, res.output
        lines = out.read_text().strip().splitlines()
        bases = lines[0].split("\t")[1:]
        rows = [dict(zip(bases, ln.split("\t")[1:])) for ln in lines[1:]]
        for row, col in zip(rows, m.scores):
            for base, score in col:
                assert float(row[base]) == pytest.approx(score)

    def test_figure_still_written_when_requested(self, tmp_path):
        pytest.importorskip("matplotlib", reason="plot extra not installed")
        mat = tmp_path / "m.tsv"
        fig = tmp_path / "logo.png"
        res = CliRunner().invoke(
            cli, ["logo", "-m", "ACGT", "-o", str(fig), "--matrix", str(mat)]
        )
        assert res.exit_code == 0, res.output
        assert fig.exists() and mat.exists()

    def test_requires_output_or_matrix(self):
        res = CliRunner().invoke(cli, ["logo", "-m", "ACGT"])
        assert res.exit_code != 0
        assert "Provide --output and/or --matrix" in res.output


class TestMetageneExportProfile:
    def test_export_profile_written(self, tmp_path):
        """--export-profile writes the profile matrix alongside -o."""
        pytest.importorskip("coralsnake.gtf", reason="metagene deps", exc_type=ImportError)
        gtf = _write_gtf(tmp_path)
        sites = tmp_path / "sites.tsv"
        sites.write_text("Chrom\tPos\tStrand\n1\t200\t+\n1\t250\t-\n1\t300\t+\n")
        prof = tmp_path / "profile.tsv"
        res = CliRunner().invoke(
            cli, ["metagene", "-i", str(sites), "-g", str(gtf), "-H",
                  "--meta-columns", "1,2,3", "--bins", "10",
                  "--export-profile", str(prof)]
        )
        assert res.exit_code == 0, res.output
        assert prof.exists()
        header = prof.read_text().splitlines()[0].split("\t")
        assert header[0] == "feature_type"
        assert "feature_midpoint" in header
        assert any(h.startswith("count") for h in header)

    def test_output_not_required_when_exporting(self, tmp_path):
        """-o/--output is no longer mandatory if --export-profile is given."""
        pytest.importorskip("coralsnake.gtf", reason="metagene deps", exc_type=ImportError)
        gtf = _write_gtf(tmp_path)
        sites = tmp_path / "sites.tsv"
        sites.write_text("Chrom\tPos\tStrand\n1\t200\t+\n")
        prof = tmp_path / "profile.tsv"
        res = CliRunner().invoke(
            cli, ["metagene", "-i", str(sites), "-g", str(gtf), "-H",
                  "--meta-columns", "1,2,3", "--bins", "5",
                  "--export-profile", str(prof)]
        )
        assert res.exit_code == 0, res.output
        assert prof.exists()

    def test_still_errors_without_any_output(self, tmp_path):
        """With no -o/-s/--export-profile the command still fails loudly."""
        gtf = _write_gtf(tmp_path)
        sites = tmp_path / "sites.tsv"
        sites.write_text("Chrom\tPos\tStrand\n1\t200\t+\n")
        res = CliRunner().invoke(
            cli, ["metagene", "-i", str(sites), "-g", str(gtf), "-H",
                  "--meta-columns", "1,2,3"]
        )
        assert res.exit_code != 0
        assert "Output file is required" in res.output
