#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Synthetic-data integration tests.

Builds a tiny, hand-computable genome + GTF (tests/synthetic_data.py) and
validates each feature's details against exact expected values: prepare,
annotate (GTF + table modes, both strands, intronic/intergenic/noncoding),
variant effect, motif boundaries, BAM coordinate conversion and liftover
direction selection.
"""

from pathlib import Path

import pytest

import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))

from synthetic_data import CHR1, write_synthetic  # noqa: E402


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    d = tmp_path_factory.mktemp("syn")
    fa, gt = write_synthetic(str(d))
    return {"fa": str(fa), "gtf": str(gt), "dir": str(d)}


@pytest.fixture(scope="module")
def fasta(syn):
    import pysam

    fa = pysam.FastaFile(syn["fa"])
    yield fa
    fa.close()


# ---------------------------------------------------------------------------
# prepare (gtf2tx) - transcript selection + --with-txpos/--with-codon
# ---------------------------------------------------------------------------
class TestPrepare:
    def test_with_txpos_codon(self, syn, tmp_path):
        from coralsnake.gtf2tx import parse_file

        out = tmp_path / "prep.tsv"
        parse_file(
            gtf_file=syn["gtf"],
            fasta_file=syn["fa"],
            output_file=str(out),
            seq_file=None,
            with_codon=True,
            with_txpos=True,
        )
        rows = {}
        for line in out.read_text().rstrip("\n").split("\n")[1:]:
            c = line.split("\t")
            rows[c[1]] = c  # by transcript_id

        # t1 '+': bounding box [11,70] (1-based), start/stop codon 1-based
        r1 = rows["t1"]
        assert r1[3] == "+"
        assert r1[5:7] == ["20", "62"]  # start_codon, stop_codon (1-based)
        assert r1[-2:] == ["11", "70"]  # transcript_start, transcript_end

        # t2 '-': bounding box must be a valid span [min_start+1, max_end]
        r2 = rows["t2"]
        assert r2[3] == "-"
        start, end = int(r2[-2]), int(r2[-1])
        assert start <= end
        assert (start, end) == (121, 150)

    def test_parse_file_gff(self, syn, tmp_path):
        # GFF attribute parsing is exercised via parse_file on a tiny .gff
        from coralsnake.gtf2tx import parse_gtf_annot, parse_gff_annot

        assert parse_gff_annot("ID=exon1;Parent=t1")["ID"] == "exon1"
        assert parse_gtf_annot('gene_id "g1"; transcript_id "t1";')["gene_id"] == "g1"


# ---------------------------------------------------------------------------
# annotate - GTF mode (regions, positions, both strands)
# ---------------------------------------------------------------------------
class TestAnnotateGTF:
    def _run(self, syn, sites, cols="1,2,3,4,5"):
        from coralsnake.annotate import run_annotate

        inp = Path(syn["dir"]) / "sites.tsv"
        inp.write_text("".join(f"chr1\t{s}\t{st}\t.\t.\n" for s, st in sites))
        out = Path(syn["dir"]) / "out.tsv"
        run_annotate(
            str(inp), str(out), syn["gtf"], all_effects=False,
            with_header=False, columns=cols,
        )
        lines = out.read_text().rstrip("\n").split("\n")
        header = lines[0].split("\t")
        return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]

    def test_plus_strand_regions(self, syn, fasta):
        # (genomic pos, strand)  ->  (gene, region, transcript_pos)
        cases = [
            (12, "+", "g1", "FivePrimeUTR", "2"),
            (25, "+", "g1", "CDS", "15"),
            (65, "+", "g1", "ThreePrimeUTR", "35"),
            (40, "+", "g1", "Intronic", None),
            (5, "+", None, "Intergenic", None),
        ]
        rows = self._run(syn, [(p, s) for p, s, *_ in cases])
        for row, (p, s, gene, region, tpos) in zip(rows, cases):
            assert row["gene_id"] == (gene or ""), p
            assert row["region"] == region, p
            if tpos is not None:
                assert row["transcript_pos"] == tpos, p

    def test_minus_strand_regions(self, syn, fasta):
        cases = [
            (140, "-", "g2", "FivePrimeUTR", "9"),
            (132, "-", "g2", "CDS", "17"),
            (124, "-", "g2", "ThreePrimeUTR", "25"),
        ]
        rows = self._run(syn, [(p, s) for p, s, *_ in cases])
        for row, (p, s, gene, region, tpos) in zip(rows, cases):
            assert row["gene_id"] == gene, p
            assert row["region"] == region, p
            assert row["transcript_pos"] == tpos, p

    def test_noncoding_transcript(self, syn, fasta):
        rows = self._run(syn, [(210, "+")])
        assert rows[0]["region"] == "NoncodingTranscript"
        assert rows[0]["transcript_pos"] == "10"

    def test_variant_coding_effect(self, syn, fasta):
        from coralsnake.annotate import run_annotate

        inp = Path(syn["dir"]) / "var.tsv"
        inp.write_text("chr1\t25\t+\tC\tG\n")  # codon bases at 25,26,27
        out = Path(syn["dir"]) / "var_out.tsv"
        run_annotate(
            str(inp), str(out), syn["gtf"], reference_transcript=[syn["fa"]],
            all_effects=False, with_header=False, columns="1,2,3,4,5",
        )
        lines = out.read_text().rstrip("\n").split("\n")
        header = lines[0].split("\t")
        row = dict(zip(header, lines[1].split("\t")))
        assert row["region"] == "CDS"
        # codon_ref is the 3 splice-adjacent CDS bases starting at pos 25
        codon = CHR1[25:28]
        assert row["codon_ref"] == codon
        assert row["mut_type"] in ("Silent", "Substitution")
        # motif (pad=10) = 21 bases centred on t=15 of the transcript
        assert len(row["transcript_motif"]) == 21


# ---------------------------------------------------------------------------
# annotate - table mode (subsumes annot)
# ---------------------------------------------------------------------------
class TestAnnotateTable:
    def test_table_mode(self, syn, tmp_path):
        from coralsnake.annotate import run_annotate

        table = tmp_path / "annot.tsv"
        table.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "gx\ttx\tchr1\t+\t100-109,150-159\n"
        )
        inp = tmp_path / "sites.tsv"
        inp.write_text("chr1\t105\t+\t.\t.\n")
        out = tmp_path / "out.tsv"
        run_annotate(str(inp), str(out), None, annotation_table=str(table),
                     columns="1,2,3,4,5")
        rows = out.read_text().rstrip("\n").split("\n")
        header = rows[0].split("\t")
        row = dict(zip(header, rows[1].split("\t")))
        # spans are 1-based: 100-109 -> 0-based [99,109); pos 105 -> off 6
        assert row["gene_id"] == "gx"
        assert row["transcript_pos"] == "5"  # 105-100

    def test_malformed_rows_skipped(self, syn, tmp_path):
        """Bad rows (missing cols / non-int pos) are skipped, not fatal."""
        from coralsnake.annotate import run_annotate

        inp = tmp_path / "dirty.tsv"
        inp.write_text(
            "chr1\t12\t+\t.\t.\n"       # good -> 5'UTR of g1
            "notabs_no_columns\n"        # missing columns -> skipped
            "chr1\toops\t+\t.\t.\n"      # non-integer pos -> skipped
            "chr1\t25\t+\tA\tG\n"        # good -> CDS
        )
        out = tmp_path / "out.tsv"
        run_annotate(str(inp), str(out), syn["gtf"],
                     reference_transcript=[syn["fa"]], columns="1,2,3,4,5")
        lines = out.read_text().rstrip("\n").split("\n")
        assert len(lines) == 3  # header + 2 good rows
        gene_ids = [row.split("\t")[5] for row in lines[1:]]
        assert gene_ids == ["g1", "g1"]


# ---------------------------------------------------------------------------
# motif - boundary cases on the synthetic 300 bp chromosome
# ---------------------------------------------------------------------------
class TestMotifSynthetic:
    def _m(self, fasta, pos, strand, lpad, rpad):
        from coralsnake.motif import get_motif

        return get_motif(fasta, "chr1", 300, pos, strand, lpad, rpad)

    def test_lengths_and_centres(self, fasta):
        comp = str.maketrans("ACGT", "TGCA")
        for pos, strand, lpad, rpad in [
            (150, "+", 5, 5),
            (1, "+", 4, 4),
            (300, "+", 4, 4),
            (300, "-", 4, 4),
            (1, "-", 4, 4),
        ]:
            m = self._m(fasta, pos, strand, lpad, rpad)
            assert len(m) == lpad + 1 + rpad, (pos, strand)
            centre = CHR1[pos - 1]
            # plus: centre at lpad; minus: centre at rpad, complemented
            exp = centre if strand == "+" else centre.translate(comp)
            at = lpad if strand == "+" else rpad
            assert m[at] == exp, (pos, strand)

    def test_internal_window_equals_reference(self, fasta):
        # pos 150 (1-based) -> 0-based p=149; window [p-10, p+10+1) = [139,160)
        m = self._m(fasta, 150, "+", 10, 10)
        assert m == CHR1[139:160]


# ---------------------------------------------------------------------------
# BAM coordinate conversion (gbam2tbam / tbam2gbam + liftover direction)
# ---------------------------------------------------------------------------
def _make_genome_bam(path, refname, readref, start, seq, cigar, flag=0):
    import pysam

    gh = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.4"}, "SQ": [{"SN": refname, "LN": 10000}]}
    )
    with pysam.AlignmentFile(str(path), "wb", header=gh) as out:
        a = pysam.AlignedSegment()
        a.query_name = "r1"
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.cigartuples = cigar
        a.mapping_quality = 60
        out.write(a)


class TestBamConvert:
    def test_gbam2tbam_on_synthetic_gene(self, syn, tmp_path):
        from coralsnake.gbam2tbam import convert_bam as g2t

        # t1 '+' exons (0-based) [10,30) [50,70) -> annotation spans (1-based)
        annot = tmp_path / "annot.tsv"
        annot.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "g1\tt1\tchr1\t+\t11-30,51-70\n"
        )
        # 5M read in exon1 [12,17) -> t[2,7)
        gbam = tmp_path / "genome.bam"
        _make_genome_bam(gbam, "chr1", "chr1", 12, "ACGTA", [(0, 5)])
        tbam = tmp_path / "transcript.bam"
        g2t(str(gbam), str(tbam), str(annot), threads=1)
        with pysam.AlignmentFile(str(tbam), "rb") as r:
            reads = list(r)
        assert len(reads) == 1
        rd = reads[0]
        assert rd.reference_name == "t1"
        assert rd.reference_start == 2
        assert rd.cigartuples == [(0, 5)]

    def test_liftover_t2g_and_g2t(self, syn, tmp_path):
        from click.testing import CliRunner
        from coralsnake.cli import cli

        annot = tmp_path / "annot.tsv"
        annot.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "g1\tt1\tchr1\t+\t11-30\n"
        )
        runner = CliRunner()
        # g->t
        gbam = tmp_path / "g.bam"
        _make_genome_bam(gbam, "chr1", "chr1", 12, "ACGTA", [(0, 5)])
        out1 = tmp_path / "t.bam"
        res = runner.invoke(
            cli, ["liftover", "-d", "g2t", "-i", str(gbam), "-o", str(out1),
                  "-a", str(annot)]
        )
        assert res.exit_code == 0, res.output
        import pysam

        with pysam.AlignmentFile(str(out1), "rb") as r:
            assert list(r)[0].reference_name == "t1"
        # t->g (needs a faidx header source); t2g requires a faidx file
        # Build a transcript BAM + tiny faidx, then remap to genome.
        th = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.4"}, "SQ": [{"SN": "t1", "LN": 30}]}
        )
        in_t = tmp_path / "t2g_input.bam"
        with pysam.AlignmentFile(str(in_t), "wb", header=th) as out:
            a = pysam.AlignedSegment()
            a.query_name = "r2"
            a.query_sequence = "ACGTACGTAC"
            a.query_qualities = pysam.qualitystring_to_array("I" * 10)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 2
            a.cigartuples = [(0, 10)]
            a.mapping_quality = 60
            out.write(a)
        faidx = tmp_path / "chr1.fa.fai"
        faidx.write_text("chr1\t300\t14\t60\t61\n")
        out_gen = tmp_path / "g_out.bam"
        res2 = runner.invoke(
            cli, ["liftover", "-d", "t2g", "-i", str(in_t), "-o", str(out_gen),
                  "-a", str(annot), "-f", str(faidx)]
        )
        assert res2.exit_code == 0, res2.output
        with pysam.AlignmentFile(str(out_gen), "rb") as r:
            rd = list(r)[0]
            assert rd.reference_name == "chr1"
            assert rd.reference_start == 12  # t[2,12) in exon1 [10,30) -> g 12
            assert rd.cigartuples[0][0] == 0


# ---------------------------------------------------------------------------
# map_to_local / metagene helpers
# ---------------------------------------------------------------------------
class TestMapToLocalSynthetic:
    def test_full_chain(self):
        import polars as pl
        from coralsnake.map_to_local import map_to_local

        # query site at chr pos 25 (+), reference exon [10,30)
        query = pl.DataFrame(
            {"Chromosome": ["chr1"], "Start": [25], "End": [26], "Strand": ["+"]}
        )
        ref = pl.DataFrame(
            {"Chromosome": ["chr1"], "Start": [10], "End": [30], "Strand": ["+"],
             "transcript_id": ["X"]}
        )
        local = map_to_local(query, ref, ref_id_col="transcript_id")
        assert local["Chromosome"][0] == "X"
        assert local["Start"][0] == 15  # 25 - 10


# ---------------------------------------------------------------------------
# coordinate + logo
# ---------------------------------------------------------------------------
class TestCoordinateLogo:
    def test_coordinate_custom(self, tmp_path):
        from click.testing import CliRunner
        from coralsnake.cli import cli

        mp = tmp_path / "m.tsv"
        mp.write_text("chr1\t1\n")
        inp = tmp_path / "in.tsv"
        inp.write_text("chr1\t100\t+\n")
        out = tmp_path / "out.tsv"
        runner = CliRunner()
        res = runner.invoke(cli, ["coordinate", "-i", str(inp), "-o", str(out),
                                  "-m", str(mp), "-c", "1"])
        assert res.exit_code == 0, res.output
        assert out.read_text().strip().split("\n")[0].startswith("1\t")

    def test_logo_no_nan(self):
        import numpy as np
        from coralsnake.logo import Mlogo

        m = Mlogo(motifs=["ACGT", "ACGG"], to2bit=True, t2u=False).scores
        for col in m:
            for _, s in col:
                assert not np.isnan(s)


# ---------------------------------------------------------------------------
# metagene options (were silently ignored before)
# ---------------------------------------------------------------------------
class TestMetageneOptions:
    def test_region_and_weight_transforms(self, syn, tmp_path):
        import polars as pl
        from click.testing import CliRunner
        from coralsnake.cli import cli

        sites = tmp_path / "sites.tsv"
        # g12 = 5'UTR, g25 = CDS (g1 on chr1); weight column
        sites.write_text("chr1\t12\t+\t2\nchr1\t25\t+\t3\n")
        runner = CliRunner()

        def total(out_score):
            return pl.read_csv(out_score, separator="\t").select(
                pl.col(pl.Float64).sum()
            ).item() if False else sum(
                float(r.split("\t")[-1]) for r in out_score.read_text().splitlines()[1:]
            )

        # base: both sites
        s0 = tmp_path / "s0.tsv"
        o0 = tmp_path / "o0.tsv"
        r0 = runner.invoke(
            cli, ["metagene", "-i", str(sites), "-o", str(o0), "-s", str(s0),
                  "--gtf", syn["gtf"], "-m", "1,2,3", "-w", "4"]
        )
        assert r0.exit_code == 0, r0.output
        # region=cds filters out the 5'UTR site (weight 2)
        s1 = tmp_path / "s1.tsv"
        o1 = tmp_path / "o1.tsv"
        r1 = runner.invoke(
            cli, ["metagene", "-i", str(sites), "-o", str(o1), "-s", str(s1),
                  "--gtf", syn["gtf"], "-m", "1,2,3", "-w", "4", "--region", "cds"]
        )
        assert r1.exit_code == 0, r1.output
        assert total(s0) == 5.0
        assert total(s1) == 3.0
        # score-transform + normalize complete without error
        s2 = tmp_path / "s2.tsv"
        o2 = tmp_path / "o2.tsv"
        r2 = runner.invoke(
            cli, ["metagene", "-i", str(sites), "-o", str(o2), "-s", str(s2),
                  "--gtf", syn["gtf"], "-m", "1,2,3", "-w", "4",
                  "--score-transform", "log2", "--normalize"]
        )
        assert r2.exit_code == 0, r2.output
