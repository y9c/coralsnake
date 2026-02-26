"""Tests for coralsnake.annot – site annotation with transcript positions."""

from pathlib import Path

import pytest

run_annot = pytest.importorskip("coralsnake.annot", reason="annot deps not available", exc_type=ImportError).run_annot


def _write(p: Path, content: str) -> str:
    p.write_text(content)
    return str(p)


# ---------------------------------------------------------------------------
# run_annot basics
# ---------------------------------------------------------------------------
class TestRunAnnot:
    def _annot(self, tmp_path):
        return _write(
            tmp_path / "annot.tsv",
            "chrom\tstrand\tspans\tgene_id\ttranscript_id\n"
            "chr1\t+\t10-20,30-40\tGENE1\tTX1\n"
            "chr1\t-\t50-60\tGENE2\tTX2\n",
        )

    def test_plus_strand_first_base(self, tmp_path):
        annot = self._annot(tmp_path)
        inp = _write(tmp_path / "sites.tsv", "chr1\t10\t+\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot, keep_na=True)
        line = Path(out).read_text().strip().split("\n")[0]
        assert line == "chr1\t10\t+\tGENE1\tTX1\t0"

    def test_plus_strand_second_exon(self, tmp_path):
        annot = self._annot(tmp_path)
        inp = _write(tmp_path / "sites.tsv", "chr1\t35\t+\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot, keep_na=True)
        line = Path(out).read_text().strip().split("\n")[0]
        assert line == "chr1\t35\t+\tGENE1\tTX1\t16"

    def test_no_hit_gives_na(self, tmp_path):
        annot = self._annot(tmp_path)
        inp = _write(tmp_path / "sites.tsv", "chr1\t5\t+\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot, keep_na=True)
        line = Path(out).read_text().strip().split("\n")[0]
        assert line == "chr1\t5\t+\t.\t.\t."

    def test_keep_na_false(self, tmp_path):
        annot = self._annot(tmp_path)
        inp = _write(tmp_path / "sites.tsv", "chr1\t5\t+\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot, keep_na=False)
        assert Path(out).read_text().strip() == ""

    def test_minus_strand_position(self, tmp_path):
        """Verify the off-by-one fix: exon_end - 1 - position + exon_shift."""
        annot = self._annot(tmp_path)
        inp = _write(tmp_path / "sites.tsv", "chr1\t60\t-\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot, keep_na=True)
        line = Path(out).read_text().strip().split("\n")[0]
        assert line == "chr1\t60\t-\tGENE2\tTX2\t0"

    def test_minus_strand_last_base(self, tmp_path):
        annot = self._annot(tmp_path)
        inp = _write(tmp_path / "sites.tsv", "chr1\t51\t-\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot, keep_na=True)
        line = Path(out).read_text().strip().split("\n")[0]
        assert line == "chr1\t51\t-\tGENE2\tTX2\t9"

    def test_collapse_annot(self, tmp_path):
        annot_path = _write(
            tmp_path / "annot.tsv",
            "chrom\tstrand\tspans\tgene_id\ttranscript_id\n"
            "chr1\t+\t10-20\tG1\tT1\n"
            "chr1\t+\t10-20\tG2\tT2\n",
        )
        inp = _write(tmp_path / "sites.tsv", "chr1\t15\t+\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot_path, keep_na=True, collapse_annot=True)
        line = Path(out).read_text().strip().split("\n")[0]
        assert "G1,G2" in line or "G2,G1" in line

    def test_add_count(self, tmp_path):
        annot_path = _write(
            tmp_path / "annot.tsv",
            "chrom\tstrand\tspans\tgene_id\ttranscript_id\n"
            "chr1\t+\t10-20\tG1\tT1\n"
            "chr1\t+\t10-20\tG2\tT2\n",
        )
        inp = _write(tmp_path / "sites.tsv", "chr1\t15\t+\n")
        out = str(tmp_path / "out.tsv")
        run_annot(inp, out, annot_path, keep_na=True, add_count=True)
        lines = Path(out).read_text().strip().split("\n")
        assert all(line.endswith("\t2") for line in lines)
