#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the unified site/variant annotation engine (coralsnake.annotate).

`annotate` merges the logic of the former ``annot`` (site labeling) and
``effect`` (variant effect) commands into a single GTF-based engine with one
fixed output schema.
"""

from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent / "data"
GTF = DATA / "R64-1-1.release57.gtf"
GENOME = DATA / "R64-1-1.fa"


class TestAnnotationSchema:
    def test_header_is_fixed(self):
        from coralsnake.annotate import Annotation

        assert Annotation.header() == [
            "gene_id",
            "transcript_id",
            "transcript_pos",
            "region",
            "gene_pos",
            "transcript_strand",
            "mut_type",
            "transcript_motif",
            "coding_pos",
            "codon_ref",
            "aa_pos",
            "aa_ref",
            "distance2splice",
        ]

    def test_annotation_columns_match_header(self):
        from coralsnake.annotate import Annotation

        a = Annotation(
            gene_id="G", transcript_id="T", transcript_pos=3, region="CDS"
        )
        assert len(a.to_list()) == len(Annotation.header())


class TestAnnotateSite:
    def test_pure_site_gets_region_without_fasta(self):
        """A bare site (no ref/alt, no FASTA) still yields region + position."""
        if not GTF.exists():
            pytest.skip("R64 GTF not present")
        from coralsnake.annotate import build_transcript_index  # noqa: F401
        from coralsnake.effect import build_transcript_index as bti, Site
        from coralsnake.annotate import _annotate_site

        by_chrom = bti(str(GTF))
        site = Site(chrom="I", pos=74019, strand="+", ref="-", alt="N")
        anns = _annotate_site(site, by_chrom)
        top = max(anns, key=lambda a: a.region or "")
        assert top.gene_id == "YAL037W"
        assert top.transcript_pos == 0
        assert top.region == "CDS"


class TestAnnotateEndToEnd:
    def test_variant_mode(self, tmp_path):
        if not (GTF.exists() and GENOME.exists()):
            pytest.skip("R64 GTF/FASTA not present")
        from coralsnake.annotate import run_annotate

        inp = tmp_path / "sites.tsv"
        inp.write_text("I\t74019\t+\tA\tG\n")
        out = tmp_path / "out.tsv"
        run_annotate(
            str(inp),
            str(out),
            str(GTF),
            reference_transcript=[str(GENOME)],
            npad=10,
            strandness=True,
            all_effects=True,
            with_header=False,
            columns="1,2,3,4,5",
        )
        lines = out.read_text().rstrip("\n").split("\n")
        assert len(lines) == 2  # header + 1 row
        header = lines[0].split("\t")
        assert "region" in header
        row = dict(zip(header, lines[1].split("\t")))
        assert row["gene_id"] == "YAL037W"
        assert row["region"] == "CDS"
        assert row["transcript_pos"] == "0"
        assert row["mut_type"] == "Substitution"
        assert row["codon_ref"] in ("ATG", "ACC")
        assert row["aa_ref"] in ("M", "N")

    def test_site_mode_without_fasta(self, tmp_path):
        if not GTF.exists():
            pytest.skip("R64 GTF not present")
        from coralsnake.annotate import run_annotate

        inp = tmp_path / "sites.tsv"
        inp.write_text("I\t74019\t+\t.\t.\n")
        out = tmp_path / "out.tsv"
        run_annotate(
            str(inp),
            str(out),
            str(GTF),
            npad=10,
            strandness=True,
            all_effects=False,
            with_header=False,
            columns="1,2,3,4,5",
        )
        lines = out.read_text().rstrip("\n").split("\n")
        header = lines[0].split("\t")
        row = dict(zip(header, lines[1].split("\t")))
        assert row["region"] == "CDS"
        assert row["gene_id"] == "YAL037W"
        # No FASTA -> coding columns stay empty.
        assert row["codon_ref"] == ""
        assert row["transcript_motif"] == ""

    def test_intergenic(self, tmp_path):
        if not GTF.exists():
            pytest.skip("R64 GTF not present")
        from coralsnake.annotate import run_annotate

        inp = tmp_path / "sites.tsv"
        inp.write_text("I\t100\t+\t.\t.\n")
        out = tmp_path / "out.tsv"
        run_annotate(
            str(inp),
            str(out),
            str(GTF),
            all_effects=False,
            with_header=False,
            columns="1,2,3,4,5",
        )
        lines = out.read_text().rstrip("\n").split("\n")
        header = lines[0].split("\t")
        row = dict(zip(header, lines[1].split("\t")))
        assert row["region"] == "Intergenic"
        assert row["mut_type"] == "Intergenic"
