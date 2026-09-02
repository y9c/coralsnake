#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Tests for the linked genome FASTA download (coralsnake.download.
# download_genome / ensure_reference_*) and the reference-name resolution in
# the CLI (annotate -g/-f, liftover -a/-f, motif -f, reference genome).
#
# A fake built-in reference "FAKE" is installed via monkeypatch: its parquet
# is pre-seeded into a temp XDG cache and its genome link points at a local
# file:// FASTA, so no network access is needed.

import gzip
import os
import time

import pytest

pytest.importorskip("polars")

import polars as pl  # noqa: E402

from coralsnake.config import BUILTIN_REFERENCES  # noqa: E402
from coralsnake.config import GENOME_URLS  # noqa: E402

GENOME_FA = (
    ">I contig one\n"
    "ACGTACGTACGTACGTACGT\n"
    "ACGTACGTACGTACGTACGT\n"
    ">II contig two\n"
    "TTGGCCAA\n"
)

FAKE_GTF = (
    "##gtf-version 2.2\n"
    "#!genome-build FAKE\n"
    'I\ttest\tgene\t10\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; '
    'gene_name "Foo"; gene_biotype "protein_coding";\n'
    'I\ttest\ttranscript\t10\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; '
    'transcript_biotype "protein_coding";\n'
    'I\ttest\texon\t10\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; '
    'exon_number "1"; gene_name "Foo";\n'
)


@pytest.fixture()
def fake_ref(tmp_path, monkeypatch):
    """Install a FAKE built-in reference backed by local files.

    Returns the cache dir (where the parquet and derived artifacts live).
    """
    from coralsnake.gtf import prepare_exon_ref

    cache = tmp_path / "cache" / "coralsnake"
    cache.mkdir(parents=True)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))

    gtf = tmp_path / "fake.gtf"
    gtf.write_text(FAKE_GTF)
    prepare_exon_ref(str(gtf)).write_parquet(cache / "FAKE.parquet")

    gz = tmp_path / "genome.fa.gz"
    gz.write_bytes(gzip.compress(GENOME_FA.encode()))

    monkeypatch.setitem(
        BUILTIN_REFERENCES,
        "FAKE",
        {
            "parquet_file": "FAKE.parquet",
            "source_file": "Fake/raw/fake.gtf.gz",
            "source_url": "https://example.invalid/fake.gtf.gz",
            "description": "Fake test reference",
        },
    )
    monkeypatch.setitem(GENOME_URLS, "FAKE", gz.as_uri())
    return cache


class TestDownloadGenome:
    def test_fetches_decompresses_and_indexes(self, fake_ref):
        from coralsnake.download import download_genome

        fa = download_genome("FAKE", silent=True)
        assert fa.name == "FAKE.fa"
        assert fa.parent == fake_ref / "genomes"
        assert fa.exists()
        assert fa.with_name(fa.name + ".fai").exists()
        assert (fake_ref / "genomes" / "FAKE.fa.gz").exists()

        import pysam

        seq = pysam.FastaFile(str(fa)).fetch("I", 0, 20)
        assert seq == "ACGTACGTACGTACGTACGT"

    def test_idempotent(self, fake_ref):
        from coralsnake.download import download_genome

        fa = download_genome("FAKE", silent=True)
        mtime = fa.stat().st_mtime_ns
        download_genome("FAKE", silent=True)
        assert fa.stat().st_mtime_ns == mtime

    def test_unknown_reference(self, fake_ref):
        from coralsnake.download import download_genome

        with pytest.raises(ValueError):
            download_genome("NO-SUCH-REF", silent=True)

    def test_header_check_reports_missing(self, fake_ref):
        from coralsnake.download import _check_genome_headers, download_genome

        fa = download_genome("FAKE", silent=True)
        # partial overlap: I is in the FASTA, NOPE is not
        pl.DataFrame({"Chromosome": ["I", "NOPE"]}).write_parquet(
            fake_ref / "FAKE.parquet"
        )
        assert _check_genome_headers("FAKE", fa) == {"NOPE"}
        # no overlap at all: hard failure
        pl.DataFrame({"Chromosome": ["NOPE2"]}).write_parquet(fake_ref / "FAKE.parquet")
        with pytest.raises(RuntimeError, match="shares no contig names"):
            _check_genome_headers("FAKE", fa)

    def test_download_rejects_mismatched_genome(self, fake_ref):
        from coralsnake.download import download_genome

        pl.DataFrame({"Chromosome": ["NOPE2"]}).write_parquet(fake_ref / "FAKE.parquet")
        with pytest.raises(RuntimeError):
            download_genome("FAKE", silent=True)


class TestEnsureDerived:
    def test_table_and_gtf_cached(self, fake_ref):
        from coralsnake.download import ensure_reference, ensure_reference_gtf
        from coralsnake.download import ensure_reference_table

        pq = ensure_reference("FAKE")
        assert pq == fake_ref / "FAKE.parquet"

        table = ensure_reference_table("FAKE")
        assert table == fake_ref / "FAKE.table.tsv"
        header = table.read_text().splitlines()[0].split("\t")
        assert header[:5] == ["gene_id", "transcript_id", "chrom", "strand", "spans"]

        gtf = ensure_reference_gtf("FAKE")
        assert gtf == fake_ref / "FAKE.gtf"
        assert "transcript" in gtf.read_text()

    def test_stale_table_reexported(self, fake_ref):
        from coralsnake.download import ensure_reference_table

        table = ensure_reference_table("FAKE")
        mtime = table.stat().st_mtime_ns
        # pretend the parquet was updated (newer than the derived table)
        new_m = time.time() + 10
        os.utime(fake_ref / "FAKE.parquet", (new_m, new_m))
        table2 = ensure_reference_table("FAKE")
        assert table2 == table
        assert table.stat().st_mtime_ns > mtime


class TestCli:
    def test_reference_genome(self, fake_ref):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        res = CliRunner().invoke(cli, ["reference", "genome", "FAKE"])
        assert res.exit_code == 0, res.output
        assert (fake_ref / "genomes" / "FAKE.fa").exists()

    def test_reference_download_with_genome(self, fake_ref):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        res = CliRunner().invoke(
            cli, ["reference", "download", "FAKE", "--with-genome"]
        )
        assert res.exit_code == 0, res.output
        assert (fake_ref / "genomes" / "FAKE.fa").exists()

    def test_reference_genome_unknown_name(self, fake_ref):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        res = CliRunner().invoke(cli, ["reference", "genome", "NO-SUCH-REF"])
        assert res.exit_code != 0

    def test_annotate_by_reference_name(self, fake_ref, tmp_path):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        sites = tmp_path / "sites.tsv"
        sites.write_text("chrom\tpos\tstrand\nI\t30\t+\n")
        out = tmp_path / "out.tsv"
        res = CliRunner().invoke(
            cli,
            ["annotate", "-i", str(sites), "-H", "-g", "FAKE", "-o", str(out)],
        )
        assert res.exit_code == 0, res.output
        rows = out.read_text().strip().splitlines()
        assert len(rows) == 2
        assert "g1" in rows[1]
        # the derived GTF is cached next to the parquet
        assert (fake_ref / "FAKE.gtf").exists()

    def test_annotate_fasta_by_reference_name(self, fake_ref, tmp_path):
        """-f <name> fetches the linked genome FASTA and runs GTF mode."""
        from click.testing import CliRunner

        from coralsnake.cli import cli

        sites = tmp_path / "sites.tsv"
        sites.write_text("chrom\tpos\tstrand\nI\t30\t+\n")
        out = tmp_path / "out.tsv"
        res = CliRunner().invoke(
            cli,
            [
                "annotate",
                "-i",
                str(sites),
                "-H",
                "-g",
                "FAKE",
                "-f",
                "FAKE",
                "-o",
                str(out),
            ],
        )
        assert res.exit_code == 0, res.output
        assert (fake_ref / "genomes" / "FAKE.fa").exists()
        assert (fake_ref / "genomes" / "FAKE.fa.fai").exists()

    def test_liftover_table_by_reference_name(self, fake_ref, tmp_path):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        sites = tmp_path / "sites.tsv"
        sites.write_text("Chrom\tPos\tStrand\nI\t30\t+\n")
        out = tmp_path / "out.tsv"
        res = CliRunner().invoke(
            cli,
            [
                "liftover",
                "--table",
                "-d",
                "g2t",
                "-a",
                "FAKE",
                "-i",
                str(sites),
                "-o",
                str(out),
            ],
        )
        assert res.exit_code == 0, res.output
        rows = out.read_text().strip().splitlines()
        assert rows[0].endswith("Gene\tGenePos")
        assert rows[1].split("\t")[-2:] == ["g1", "21"]
        assert (fake_ref / "FAKE.table.tsv").exists()

    def test_unknown_path_errors(self, fake_ref, tmp_path):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        sites = tmp_path / "sites.tsv"
        sites.write_text("chrom\tpos\tstrand\nI\t30\t+\n")
        out = tmp_path / "out.tsv"
        res = CliRunner().invoke(
            cli,
            [
                "annotate",
                "-i",
                str(sites),
                "-H",
                "-g",
                "/nonexistent/gene.gtf",
                "-o",
                str(out),
            ],
        )
        assert res.exit_code != 0
        assert "not an existing file" in res.output

    def test_explicit_file_still_works(self, fake_ref, tmp_path):
        """A real file path takes precedence over reference-name resolution."""
        from click.testing import CliRunner

        from coralsnake.cli import cli
        from coralsnake.gtf import prepare_exon_ref
        from coralsnake.ref_export import export_gtf

        gtf = fake_ref / "FAKE.gtf"
        export_gtf(prepare_exon_ref(str(tmp_path / "fake.gtf")), str(gtf))
        sites = tmp_path / "sites.tsv"
        sites.write_text("chrom\tpos\tstrand\nI\t30\t+\n")
        out = tmp_path / "out.tsv"
        res = CliRunner().invoke(
            cli,
            ["annotate", "-i", str(sites), "-H", "-g", str(gtf), "-o", str(out)],
        )
        assert res.exit_code == 0, res.output
        assert "g1" in out.read_text()
