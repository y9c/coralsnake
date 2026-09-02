"""Regression tests for the 2026-09 code audit.

Each test guards one verified bug fix. Where the audit pinned down that a
behavior was *already correct* (e.g. the SAM-convention minus-strand handling
in t2g), a ground-truth test locks the correct behavior in.
"""
import time
from pathlib import Path

import pysam
import pytest

DATA = Path(__file__).resolve().parent / "data"
GTF = DATA / "R64-1-1.release57.gtf"
GENOME = DATA / "R64-1-1.fa"


def RC(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _align(ref_id, start, seq, cigar, flag=0, mate_ref=None, mate_start=None):
    a = pysam.AlignedSegment()
    a.query_name = "r1"
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.flag = flag
    a.reference_id = ref_id
    a.reference_start = start
    a.cigartuples = cigar
    a.mapping_quality = 60
    if mate_ref is not None:
        a.next_reference_id = mate_ref
        a.next_reference_start = mate_start
    return a


# ---------------------------------------------------------------------------
# tbam2gbam (t2g)
# ---------------------------------------------------------------------------
class TestT2GRegressions:
    def _tx_bam(self, tmp_path, seq, cigar, start=0, flag=0,
                mate_ref=None, mate_start=None, ref_len=100):
        th = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.4"}, "SQ": [{"SN": "t1", "LN": ref_len}]}
        )
        bam = tmp_path / "in.bam"
        with pysam.AlignmentFile(str(bam), "wb", header=th) as out:
            out.write(_align(0, start, seq, cigar, flag,
                             mate_ref=mate_ref, mate_start=mate_start))
        annot = tmp_path / "a.tsv"
        annot.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "g1\tt1\tchr1\t+\t101-150,201-250\n"
        )
        fai = tmp_path / "g.fa.fai"
        fai.write_text("chr1\t1000\t6\t60\t61\n")
        return str(bam), str(annot), str(fai)

    def test_eq_x_ops_consume_reference(self, tmp_path):
        """CIGAR =/X must consume the reference: the intron N is inserted."""
        from coralsnake.tbam2gbam import convert_bam

        bam, annot, fai = self._tx_bam(tmp_path, "ACGTACGTAC",
                                       [(7, 5), (8, 3), (0, 2)], start=45)
        out = tmp_path / "out.bam"
        convert_bam(bam, str(out), annot, fai, threads=1)
        with pysam.AlignmentFile(str(out), "rb") as f:
            r = list(f)[0]
        # t45-49 -> g145-149 (5=), intron, t50-54 -> g200-204 (3X 2M)
        assert r.reference_start == 145
        assert r.cigarstring == "5=50N3X2M"

    def test_mate_position_remapped(self, tmp_path):
        """The mate's transcript coordinate is remapped onto the genome."""
        from coralsnake.tbam2gbam import convert_bam

        bam, annot, fai = self._tx_bam(
            tmp_path, "ACGTACGTAC", [(0, 10)], start=20, flag=0x1 | 0x2,
            mate_ref=0, mate_start=45,
        )
        out = tmp_path / "out.bam"
        convert_bam(bam, str(out), annot, fai, threads=1)
        with pysam.AlignmentFile(str(out), "rb") as f:
            r = list(f)[0]
        assert r.next_reference_name == "chr1"
        assert r.next_reference_start == 145  # t45 on '+' -> g145 (exact)

    def test_out_of_range_read_demoted_not_crash(self, tmp_path):
        """A read whose position exceeds the reference becomes unmapped."""
        from coralsnake.tbam2gbam import (parse_alignment, remap_to_genome,  # noqa: F401
                                          transcript_to_genome)
        from coralsnake.utils import Transcript, Span

        # boundary check: pos == length is out of range (was IndexError)
        tx = Transcript(gene_id="g", transcript_id="t", chrom="c", strand="+",
                        exons={1: Span(1000, 1100), 2: Span(2000, 2050)})
        with pytest.raises(ValueError):
            transcript_to_genome(tx.length, tx)

        # parse_alignment demotes a malformed read (position at the reference
        # end) to unmapped instead of raising
        gh = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.4"}, "SQ": [{"SN": "c", "LN": 10000}]}
        )
        th = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.4"}, "SQ": [{"SN": "t", "LN": 150}]}
        )
        align = pysam.AlignedSegment.fromstring("r1\t0\tt\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII", th)
        align.reference_start = tx.length  # past the end -> out-of-range remap
        annot = {"t": tx}
        out = parse_alignment(align, annot, gh)
        assert out.reference_id == -1  # demoted to unmapped


# ---------------------------------------------------------------------------
# gbam2tbam (g2t)
# ---------------------------------------------------------------------------
class TestG2TRegressions:
    def _gbam(self, tmp_path, chrom="chr1", start=137, seq=None, cigar=None,
              flag=0x10, annot_strand="-"):
        seq = seq or "ACGTA"
        cigar = cigar or [(0, 5)]
        gh = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.4"}, "SQ": [{"SN": "chr1", "LN": 10000},
                                          {"SN": "chrX", "LN": 10000}]}
        )
        bam = tmp_path / "g.bam"
        ref_id = 0 if chrom == "chr1" else 1
        with pysam.AlignmentFile(str(bam), "wb", header=gh) as out:
            out.write(_align(ref_id, start, seq, cigar, flag))
        annot = tmp_path / "a.tsv"
        annot.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            f"g1\tt1\tchr1\t{annot_strand}\t121-150\n"
        )
        return str(bam), str(annot)

    def test_minus_strand_seq_is_reverse_complemented(self, tmp_path):
        """g2t must RC the stored sequence for '-' transcripts (mirror of t2g).

        Ground truth from the SAM convention: SEQ is stored reference-forward;
        a '-' transcript reference is RC of the genomic forward strand.
        """
        from coralsnake.gbam2tbam import convert_bam

        # genome bases at [137,142) = 'ACGTA'; the molecule (mRNA) = RC of that
        genome_piece = "ACGTA"
        stored_genome = RC(genome_piece)  # what an aligner stores for flag 16
        bam, annot = self._gbam(tmp_path, seq=stored_genome)
        out = tmp_path / "t.bam"
        convert_bam(bam, str(out), annot, threads=1)
        with pysam.AlignmentFile(str(out), "rb") as f:
            r = list(f)[0]
        assert r.flag & 16 == 0
        assert r.reference_start == 8  # g141..137 -> t8..12 on '-' exon [120,150)
        assert r.query_sequence == genome_piece  # RC'd back to the mRNA orientation

    def test_cross_chromosome_read_is_skipped(self, tmp_path):
        """A read on chrX must never be assigned to a chr1 transcript."""
        from coralsnake.gbam2tbam import convert_bam

        bam, annot = self._gbam(tmp_path, chrom="chrX", start=125, flag=0)
        out = tmp_path / "t.bam"
        convert_bam(bam, str(out), annot, threads=1)
        with pysam.AlignmentFile(str(out), "rb") as f:
            assert list(f) == []


# ---------------------------------------------------------------------------
# annotate (unified site/variant annotation)
# ---------------------------------------------------------------------------
class TestAnnotateRegressions:
    def _synthetic(self, tmp_path):
        import sys

        sys.path.insert(0, str(Path(__file__).parent))
        from synthetic_data import write_synthetic

        return write_synthetic(str(tmp_path))

    def test_stop_codon_bases_are_cds(self, tmp_path):
        """The 3 bases of the stop codon are CDS, not ThreePrimeUTR."""
        from coralsnake.annotate import run_annotate

        fa, gtf = self._synthetic(tmp_path)
        # synthetic g1: stop_codon at 1-based 62-64 (t=31..33); pos 65 is 3'UTR
        inp = tmp_path / "v.tsv"
        inp.write_text("chr1\t64\t+\tT\tA\n")  # 3rd stop-codon base
        out = tmp_path / "o.tsv"
        run_annotate(str(inp), str(out), str(gtf), reference_transcript=[fa],
                     columns="1,2,3,4,5")
        row = out.read_text().rstrip("\n").split("\n")[1].split("\t")
        header = out.read_text().rstrip("\n").split("\n")[0].split("\t")
        d = dict(zip(header, row))
        assert d["region"] == "CDS"

    def test_metagene_stop_codon_bases_are_cds(self):
        """Metagene feature_type keeps the 3 stop-codon bases as CDS (matches
        test_stop_codon_bases_are_cds / effect._classify_exonic). CDS runs
        [start_codon_pos, stop_codon_pos + 3), so offsets stop+1 / stop+2 are
        CDS and the 3'UTR begins at stop+3."""
        import polars as pl
        from coralsnake.annotation import normalize_positions

        def _site(pos):
            return {
                "transcript_id": "t1",
                "transcript_length": 100,
                "start_codon_pos": 30,
                "stop_codon_pos": 60,  # first base of the stop codon -> CDS [30, 63)
                "transcript_start": pos,
                "transcript_end": pos + 1,
                "record_id": f"r{pos}",
            }

        df = pl.DataFrame([_site(61), _site(62), _site(63), _site(29)])
        _, stats, _ = normalize_positions(df, split_strategy="median", bin_number=100)
        assert stats.get("CDS") == 2          # 61, 62 = 2nd/3rd stop-codon bases
        assert stats.get("3UTR") == 1         # 63 = first base after the stop codon
        assert stats.get("5UTR") == 1         # 29 = just before the start codon

    def test_metagene_noncoding_excluded(self):
        """Sites on a noncoding transcript (no start/stop codon) are 'None',
        not CDS, and never contribute to the metagene profile bins."""
        import polars as pl
        from coralsnake.annotation import normalize_positions

        def _site(pos, start, stop):
            return {
                "transcript_id": "t1",
                "transcript_length": 100,
                "start_codon_pos": start,
                "stop_codon_pos": stop,
                "transcript_start": pos,
                "transcript_end": pos + 1,
                "record_id": f"r{pos}",
            }

        df = pl.DataFrame([_site(45, 30, 60), _site(70, None, None)])
        gene_bins, stats, _ = normalize_positions(
            df, split_strategy="median", bin_number=100
        )
        assert stats == {"CDS": 1.0, "None": 1.0}
        # only the coding site lands in a bin (total == 1.0, not 2.0)
        assert gene_bins["count"].sum() == 1.0

    def test_placeholder_ref_alt_not_a_variant(self, tmp_path):
        """'.' ref/alt placeholders must not fabricate a CDS effect."""
        from coralsnake.annotate import run_annotate

        fa, gtf = self._synthetic(tmp_path)
        inp = tmp_path / "v.tsv"
        inp.write_text("chr1\t26\t+\t.\t.\n")
        out = tmp_path / "o.tsv"
        run_annotate(str(inp), str(out), str(gtf), reference_transcript=[fa],
                     columns="1,2,3,4,5")
        lines = out.read_text().rstrip("\n").split("\n")
        d = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert d["region"] == "CDS"
        assert d["mut_type"] == "CDS"  # no ref/alt -> region, not a fabricated effect

    def test_mnp_is_complex_substitution(self, tmp_path):
        from coralsnake.annotate import run_annotate

        fa, gtf = self._synthetic(tmp_path)
        inp = tmp_path / "v.tsv"
        inp.write_text("chr1\t26\t+\tCC\tGA\n")  # 2-base substitution
        out = tmp_path / "o.tsv"
        run_annotate(str(inp), str(out), str(gtf), reference_transcript=[fa],
                     columns="1,2,3,4,5")
        lines = out.read_text().rstrip("\n").split("\n")
        d = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert d["mut_type"] == "ComplexSubstitution"

    def test_table_mode_uniform_field_count(self, tmp_path):
        """Table mode: header and every row carry the same number of fields,
        the first data row is not leaked into the header, and malformed rows
        are skipped."""
        from coralsnake.annotate import run_annotate

        table = tmp_path / "annot.tsv"
        table.write_text(
            "chrom\tstrand\tspans\tgene_id\ttranscript_id\nchr1\t+\t10-20\tG1\tT1\n"
        )
        inp = tmp_path / "s.tsv"
        inp.write_text("chr1\t15\t+\nchr1\t99\t+\nchr1\tbad\t+\n")
        out = tmp_path / "o.tsv"
        run_annotate(str(inp), str(out), None, annotation_table=str(table),
                     columns="1,2,3")
        lines = out.read_text().rstrip("\n").split("\n")
        widths = {len(row.split("\t")) for row in lines}
        assert widths == {16}  # 3 input + 13 unified columns
        assert lines[0].startswith("chrom\tpos\tstrand\tgene_id\t")  # no data leak
        assert len(lines) == 3  # header + 2 valid rows (bad row skipped)


# ---------------------------------------------------------------------------
# gtf2tx / motif / io / annot-cache
# ---------------------------------------------------------------------------
class TestGtf2TxCodons:
    def test_ensembl_style_codon_lines(self, tmp_path):
        """start/stop codon lines without exon_number are parsed."""
        from coralsnake.gtf2tx import parse_file

        gtf = tmp_path / "a.gtf"
        gtf.write_text(
            'I\tens\ttranscript\t100\t250\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'I\tens\texon\t100\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";\n'
            'I\tens\texon\t200\t250\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "2";\n'
            'I\tens\tstart_codon\t105\t107\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
            'I\tens\tstop_codon\t245\t247\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        )
        out = tmp_path / "o.tsv"
        parse_file(str(gtf), None, str(out), with_codon=True)
        row = out.read_text().rstrip("\n").split("\n")[1].split("\t")
        header = out.read_text().rstrip("\n").split("\n")[0].split("\t")
        d = dict(zip(header, row))
        assert d["start_codon"] == "105"
        assert d["stop_codon"] == "245"


class TestMotifEdge:
    def test_site_past_contig_end_is_all_n(self, tmp_path):
        from coralsnake.motif import run_motif

        fa = tmp_path / "g.fa"
        fa.write_text(">c\n" + "ACGT" * 25 + "\n")  # 100 bp
        pysam.faidx(str(fa))
        inp = tmp_path / "s.tsv"
        inp.write_text("c\t1000\t+\n")
        out = tmp_path / "o.tsv"
        run_motif(str(inp), str(out), str(fa), 3, 3, False, "1,2,3",
                  to_upper=True, wrap_site=False)
        motif = out.read_text().rstrip("\n").split("\n")[0].split("\t")[3]
        assert motif == "NNNNNNN"  # 3 + 1 + 3, all out of bounds


class TestLoadSitesHeaderNames:
    def test_canonical_header_names(self, tmp_path):
        """A header named Chromosome/Start/End/Strand must not crash."""
        from coralsnake.io import load_sites

        f = tmp_path / "s.tsv"
        f.write_text("Chromosome\tStart\tEnd\tStrand\tscore\nchr1\t0\t1\t+\t5\n")
        df = load_sites(str(f), with_header=True, meta_col_index=[0, 1, 2, 3])
        assert "Chromosome" in df.columns and "_original_Chromosome" in df.columns

    def test_out_of_range_meta_columns_friendly_error(self, tmp_path):
        from coralsnake.io import load_sites

        f = tmp_path / "s.tsv"
        f.write_text("chr1\t0\t1\t+\n")
        with pytest.raises(ValueError, match="meta columns"):
            load_sites(str(f), meta_col_index=[0, 1, 2, 9])


class TestAnnotCache:
    def test_stale_pickle_cache_is_reparsed(self, tmp_path):
        """A pickle cache older than the table must not be reused."""
        from coralsnake.annot import parse_annot_file

        tab = tmp_path / "t.tsv"
        tab.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "G1\tT1\tchr1\t+\t101-200\n"
        )
        tree, info = parse_annot_file(str(tab), cache=True)
        assert len(info) == 1
        assert (tmp_path / "t.tsv.pickle").exists()

        # regenerate the table with an extra gene, older than nothing
        time.sleep(0.01)
        tab.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "G1\tT1\tchr1\t+\t101-200\n"
            "G2\tT2\tchr1\t+\t301-400\n"
        )
        import os

        os.utime(tab, None)  # bump mtime
        tree2, info2 = parse_annot_file(str(tab), cache=True)
        assert len(info2) == 2  # fresh table parsed, stale cache ignored
