#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Round-trip tests for coralsnake.ref_export (reference parquet ->
# `prepare` table / GTF) and the download-request resolution.
#
# The reference is built from the committed test GTF, which contains real
# start/stop codon features, so the codon coordinate round-trip is exercised
# on both strands.

import pytest

pytest.importorskip("polars")

import polars as pl  # noqa: E402

from coralsnake.download import _resolve_references  # noqa: E402
from coralsnake.gtf import prepare_exon_ref  # noqa: E402
from coralsnake.ref_export import export_gtf, export_table  # noqa: E402
from coralsnake.utils import load_annotation  # noqa: E402

TEST_GTF = "R64-1-1.release57.gtf"

# Columns that must survive a parquet -> GTF -> parquet round-trip exactly.
ROUNDTRIP_COLS = [
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "gene_id",
    "transcript_id",
    "exon_number",
    "transcript_length",
    "Start_exon",
    "End_exon",
    "start_codon_pos",
    "stop_codon_pos",
    "gene_name",
    "gene_biotype",
]


@pytest.fixture
def ref(data_dir):
    """Exon reference built from the committed test GTF (v2 schema)."""
    return prepare_exon_ref(str(data_dir / TEST_GTF))


class TestExportTable:
    def test_header_matches_prepare(self, ref, tmp_path):
        out = tmp_path / "table.tsv"
        export_table(ref, str(out))
        header = out.read_text().splitlines()[0].split("\t")
        assert header[:5] == ["gene_id", "transcript_id", "chrom", "strand", "spans"]
        assert "start_codon" in header and "stop_codon" in header
        assert "transcript_start" in header and "transcript_end" in header
        # v2 identity columns are carried through
        assert "gene_name" in header and "transcript_biotype" in header

    def test_spans_match_reference(self, ref, tmp_path):
        """Per-transcript exon spans must equal the reference exons
        (1-based inclusive), regardless of in-file order."""
        out = tmp_path / "table.tsv"
        export_table(ref, str(out))
        annot = load_annotation(str(out), with_header=True)

        expected = {}
        for r in ref.select(["transcript_id", "Start", "End"]).to_dicts():
            expected.setdefault(r["transcript_id"], []).append(
                (r["Start"] + 1, r["End"])
            )
        expected = {k: sorted(v) for k, v in expected.items()}

        seen = 0
        for _gene, txs in annot.items():
            for tx_id, t in txs.items():
                got = sorted((e.start + 1, e.end) for e in t.exons.values())
                assert got == [tuple(v) for v in expected[tx_id]], tx_id
                seen += 1
        assert seen == ref["transcript_id"].n_unique()

    def test_minus_strand_spans_in_tx_order(self, ref, tmp_path):
        """Minus-strand transcripts with >=2 exons: spans must be in 5'->3'
        order, i.e. strictly decreasing genomic starts."""
        out = tmp_path / "table.tsv"
        export_table(ref, str(out))
        header = out.read_text().splitlines()[0].split("\t")
        i_tx, i_strand, i_spans = (
            header.index("transcript_id"),
            header.index("strand"),
            header.index("spans"),
        )
        multi_minus = set(
            ref.filter(
                (pl.col("Strand") == "-") & (pl.len().over("transcript_id") >= 2)
            )["transcript_id"]
            .unique()
            .to_list()
        )
        checked = 0
        for line in out.read_text().splitlines()[1:]:
            f = line.split("\t")
            if f[i_tx] not in multi_minus:
                continue
            starts = [int(s.split("-")[0]) for s in f[i_spans].split(",")]
            assert f[i_strand] == "-"
            assert starts == sorted(starts, reverse=True), f[i_tx]
            checked += 1
        assert checked > 0, "test GTF should contain multi-exon minus transcripts"

    def test_codon_columns_present(self, ref, tmp_path):
        out = tmp_path / "table.tsv"
        export_table(ref, str(out))
        header = out.read_text().splitlines()[0].split("\t")
        i_sc, i_sp = header.index("start_codon"), header.index("stop_codon")
        rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
        n_sc = sum(1 for r in rows if r[i_sc])
        n_sp = sum(1 for r in rows if r[i_sp])
        # one codon cell per transcript that has a codon position
        # (codon_pos is broadcast across a transcript's exon rows, so count
        # unique transcripts, not rows)
        n_tx_sc = (
            ref.select(["transcript_id", "start_codon_pos"])
            .unique()
            .filter(pl.col("start_codon_pos").is_not_null())
            .height
        )
        n_tx_sp = (
            ref.select(["transcript_id", "stop_codon_pos"])
            .unique()
            .filter(pl.col("stop_codon_pos").is_not_null())
            .height
        )
        assert n_sc == n_tx_sc
        assert n_sp == n_tx_sp


class TestExportGtf:
    def test_roundtrip(self, ref, tmp_path):
        """parquet -> GTF -> prepare_exon_ref must reproduce every stable
        column, including codon positions re-derived from the emitted
        genomic coordinates."""
        out_gtf = tmp_path / "ref.gtf"
        export_gtf(ref, str(out_gtf))
        ref2 = prepare_exon_ref(str(out_gtf))

        key = ["Chromosome", "Start", "End", "Strand", "gene_id", "transcript_id"]
        a = ref.select(key + ROUNDTRIP_COLS[len(key) :]).sort(key + ["Start_exon"])
        b = ref2.select(key + ROUNDTRIP_COLS[len(key) :]).sort(key + ["Start_exon"])
        assert a.equals(b)

    def test_codon_rows_present(self, ref, tmp_path):
        out_gtf = tmp_path / "ref.gtf"
        export_gtf(ref, str(out_gtf))
        feats = [
            line.split("\t")[2] for line in out_gtf.read_text().splitlines() if line
        ]
        n_sc = feats.count("start_codon")
        n_sp = feats.count("stop_codon")
        # one codon row per transcript that has a codon position (codon_pos
        # is broadcast across a transcript's exon rows, so count unique
        # transcripts, not rows)
        n_tx_sc = (
            ref.select(["transcript_id", "start_codon_pos"])
            .unique()
            .filter(pl.col("start_codon_pos").is_not_null())
            .height
        )
        n_tx_sp = (
            ref.select(["transcript_id", "stop_codon_pos"])
            .unique()
            .filter(pl.col("stop_codon_pos").is_not_null())
            .height
        )
        assert n_sc == n_tx_sc
        assert n_sp == n_tx_sp

    def test_gtf_format_sane(self, ref, tmp_path):
        out_gtf = tmp_path / "ref.gtf"
        export_gtf(ref, str(out_gtf))
        lines = out_gtf.read_text().splitlines()
        assert lines, "empty GTF"
        for line in lines:
            cols = line.split("\t")
            assert len(cols) == 9, line
            s, e = int(cols[3]), int(cols[4])
            assert 1 <= s <= e, line
            assert cols[5] == "."
            assert cols[6] in ("+", "-")
            assert cols[8].endswith(";")
        # every transcript appears on a transcript row
        tx_rows = {
            line.split("\t")[8].split('transcript_id "')[1].split('"')[0]
            for line in lines
            if line.split("\t")[2] == "transcript"
        }
        assert tx_rows == set(ref["transcript_id"].unique().to_list())


class TestCliExport:
    @pytest.fixture(autouse=True)
    def _cache_ref(self, ref, tmp_path, monkeypatch):
        from click.testing import CliRunner  # noqa: F401  (ensures import works)

        cache = tmp_path / "cache" / "coralsnake"
        cache.mkdir(parents=True)
        ref.write_parquet(str(cache / "R64-1-1.parquet"))
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
        self.tmp_path = tmp_path

    def test_reference_export_table(self, ref):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        out = self.tmp_path / "table.tsv"
        res = CliRunner().invoke(
            cli, ["reference", "export", "R64-1-1", "--table", str(out)]
        )
        assert res.exit_code == 0, res.output
        assert out.exists()
        assert "spans" in out.read_text().splitlines()[0]

    def test_reference_export_gtf(self, ref):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        out = self.tmp_path / "ref.gtf"
        res = CliRunner().invoke(
            cli, ["reference", "export", "R64-1-1", "--gtf", str(out)]
        )
        assert res.exit_code == 0, res.output
        assert out.exists()
        assert "transcript" in out.read_text()

    def test_reference_export_requires_flag(self):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        res = CliRunner().invoke(cli, ["reference", "export", "GRCh38"])
        assert res.exit_code != 0

    def test_reference_list(self):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        res = CliRunner().invoke(cli, ["reference", "list"])
        assert res.exit_code == 0, res.output
        assert "GRCh38" in res.output

    def test_metagene_export_table_deprecated(self, ref):
        """The old metagene flags still work (deprecated aliases)."""
        from click.testing import CliRunner

        from coralsnake.cli import cli

        out = self.tmp_path / "table.tsv"
        res = CliRunner().invoke(
            cli, ["metagene", "-r", "R64-1-1", "--export-table", str(out)]
        )
        assert res.exit_code == 0, res.output
        assert out.exists()
        assert "deprecated" in res.output.lower()


class TestDownloadGroups:
    def test_resolve(self):
        from coralsnake.config import BUILTIN_REFERENCES, REFERENCE_GROUPS

        assert _resolve_references("all") == list(BUILTIN_REFERENCES.keys())
        assert _resolve_references("GRCh38") == ["GRCh38"]
        assert _resolve_references("human") == REFERENCE_GROUPS["human"]
        assert _resolve_references("mouse") == REFERENCE_GROUPS["mouse"]
        with pytest.raises(ValueError):
            _resolve_references("no-such-ref")
