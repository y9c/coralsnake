#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the fused variant commands (coralsnake.motif / .coordinate / .effect)."""

from pathlib import Path

import pytest

from coralsnake.effect import Annot, Site, expand_base, reverse_base

DATA = Path(__file__).resolve().parent / "data"
# R64-1-1.fa is bundled with the coralsnake tests (has contig "I" and "II").
GTF = DATA / "R64-1-1.release57.gtf"
GENOME = DATA / "R64-1-1.fa"


class TestConstants:
    def test_expand_base(self):
        assert expand_base("N") == ["G", "A", "T", "C"]
        assert expand_base("A") == ["A"]
        assert expand_base("R") == ["A", "G"]

    def test_reverse_base(self):
        assert reverse_base("AAGC") == "GCTT"
        assert reverse_base("") == ""

    def test_annot_column_order(self):
        # Output column order is fixed / identical to the standalone package.
        assert Annot().get_names() == [
            "mut_type",
            "gene_type",
            "gene_name",
            "gene_pos",
            "transcript_name",
            "transcript_pos",
            "transcript_motif",
            "transcript_strand",
            "coding_pos",
            "codon_ref",
            "aa_pos",
            "aa_ref",
            "distance2splice",
        ]


class TestMotif:
    @pytest.fixture
    def fasta(self):
        if not GENOME.exists():
            pytest.skip("R64 FASTA not present")
        import pysam

        fa = pysam.FastaFile(str(GENOME))
        yield fa
        fa.close()

    def test_plus_strand(self, fasta):
        chrom = fasta.references[0]
        from coralsnake.motif import get_motif

        seq = get_motif(fasta, chrom, fasta.get_reference_length(chrom), 100, "+", 3, 3)
        assert len(seq) == 7
        assert set(seq.upper()) <= set("ACGTN")  # only nucleotide bases

    def test_minus_strand_is_reverse(self, fasta):
        chrom = fasta.references[0]
        from coralsnake.motif import get_motif

        plus = get_motif(
            fasta, chrom, fasta.get_reference_length(chrom), 100, "+", 3, 3
        )
        minus = get_motif(
            fasta, chrom, fasta.get_reference_length(chrom), 100, "-", 3, 3
        )
        # minus should be the reverse-complement of plus (with N padding).
        from coralsnake.utils import reverse_complement

        assert reverse_complement(plus) == minus

    def test_padding_out_of_bounds(self, fasta):
        chrom = fasta.references[0]
        from coralsnake.motif import get_motif

        seq = get_motif(fasta, chrom, fasta.get_reference_length(chrom), 1, "+", 3, 3)
        assert len(seq) == 7
        assert seq[:3] == "NNN"  # left-padded past the chromosome start


class TestCoordinate:
    @pytest.fixture
    def chrom_map(self, tmp_path):
        p = tmp_path / "chrom_map.tsv"
        p.write_text("chr1\t1\nchrX\tX\nchrM\tMT\n")
        return str(p)

    def test_custom_mapping(self, tmp_path, chrom_map):
        from coralsnake.coordinate import run_coordinate

        inp = tmp_path / "in.tsv"
        inp.write_text("chr1\t100\t+\nchrM\t50\t-\n")
        out = tmp_path / "out.tsv"
        run_coordinate(str(inp), str(out), chrom_map, None, "1", False, False)
        lines = out.read_text().strip().split("\n")
        assert lines[0] == "1\t100\t+"
        assert lines[1] == "MT\t50\t-"

    def test_builtin_u2e(self, tmp_path):
        from coralsnake.coordinate import run_coordinate

        inp = tmp_path / "in.tsv"
        inp.write_text("chr1\t100\nchrM\t50\n")
        out = tmp_path / "out.tsv"
        run_coordinate(str(inp), str(out), None, "U2E", "1", False, False)
        lines = out.read_text().strip().split("\n")
        assert lines[0] == "1\t100"
        assert lines[1] == "MT\t50"


class TestEffect:
    def test_refine_cds_effect_position_aware(self):
        # Mutation is applied at the exact base, not the first occurrence.
        from coralsnake.effect import _refine_cds_effect

        # codon 'AAA', mutate the 2nd base (idx 1) to G → 'AGA' (Arg)
        assert _refine_cds_effect("AAA", "A", "G", 1) == "Substitution"
        # mutate the 1st base (idx 0) to G → 'GAA' (Glu)
        assert _refine_cds_effect("AAA", "A", "G", 0) == "Substitution"
        # no change → silent
        assert _refine_cds_effect("AAA", "A", "A", 0) == "Silent"
        # synonymous (different codon, same amino acid) → silent
        assert _refine_cds_effect("TTA", "A", "G", 2) == "Silent"  # TTA→TTG, both Leu
        assert _refine_cds_effect("GGA", "A", "C", 2) == "Silent"  # GGA→GGC, both Gly
        # introducing a DNA stop codon → premature stop (GAG→TAG)
        assert _refine_cds_effect("GAG", "G", "T", 0) == "PrematureStop"
        # net length change not a multiple of three → frameshift
        assert _refine_cds_effect("ATG", "AT", "A", 0) == "FrameShift"
        # net length change a multiple of three → in-frame indel
        assert _refine_cds_effect("ATG", "A", "AAAA", 0) == "InFrameIndel"

    def test_runs_end_to_end(self, tmp_path):
        if not GTF.exists():
            pytest.skip("R64 GTF not present")
        from coralsnake.effect import run_effect

        inp = tmp_path / "sites.tsv"
        inp.write_text("I\t74019\t+\tA\tG\n")
        out = tmp_path / "out.tsv"
        run_effect(
            str(inp),
            str(out),
            reference_gtf=str(GTF),
            reference_transcript=[str(GENOME)],
            reference_protein=[],
            npad=10,
            strandness=True,
            all_effects=True,
            pU_mode=False,
            with_header=False,
            columns="1,2,3,4,5",
        )
        lines = out.read_text().strip().split("\n")
        # header + 1 row
        assert len(lines) == 2
        assert lines[0].startswith("chrom\tpos")
        row = lines[1].split("\t")
        assert row[7] == "YAL037W"  # gene_name
        assert row[14] in ("ATG", "ACC")  # codon_ref inside CDS


class TestSite:
    def test_to_list(self):
        site = Site(chrom="chr1", pos=100, strand="+", ref="A", alt="G")
        assert site.to_list() == ["chr1", 100, "+", "A", "G"]
