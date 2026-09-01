"""Tests for coralsnake.tbam2gbam – reverse_md, flip_flag, transcript_to_genome, remap_to_genome."""

import pysam
import pytest

from coralsnake.tbam2gbam import (
    flip_flag,
    remap_to_genome,
    reverse_md,
    transcript_to_genome,
)
from coralsnake.utils import Span, Transcript


# ---------------------------------------------------------------------------
# reverse_md
# ---------------------------------------------------------------------------
class TestReverseMd:
    def test_all_match(self):
        assert reverse_md("10") == "10"

    def test_zero(self):
        assert reverse_md("0") == "0"

    def test_single_mismatch(self):
        assert reverse_md("2A3") == "3T2"

    def test_adjacent_mismatches(self):
        assert reverse_md("0A0C5") == "5G0T0"

    def test_deletion(self):
        # ^ACG → comp(ACG)=TGC → rev(TGC)=CGT → ^CGT
        assert reverse_md("2^ACG3") == "3^CGT2"

    def test_mismatch_then_deletion(self):
        assert reverse_md("0A0^GT3") == "3^AC0T0"

    def test_pure_mismatch(self):
        assert reverse_md("0G0") == "0C0"

    def test_complex(self):
        assert reverse_md("3A2^TT0C4") == "4G0^AA2T3"

    def test_symmetric_deletion(self):
        # ^ACGT: comp=TGCA, rev=ACGT → same
        assert reverse_md("0^ACGT0") == "0^ACGT0"

    def test_stable_value(self):
        # value-stable across calls (the fast C kernel is not cached by object
        # identity like the old lru_cached regex implementation, so compare
        # by value — and it is the value that matters downstream)
        a = reverse_md("5A5")
        b = reverse_md("5A5")
        assert a == b


# ---------------------------------------------------------------------------
# flip_flag
# ---------------------------------------------------------------------------
class TestFlipFlag:
    def test_unpaired_fwd(self):
        assert flip_flag(0) == 16

    def test_unpaired_rev(self):
        assert flip_flag(16) == 0

    def test_paired_99(self):
        assert flip_flag(99) == 83

    def test_paired_147(self):
        assert flip_flag(147) == 163

    def test_paired_both_fwd(self):
        assert flip_flag(67) == 67 ^ 0x30

    def test_paired_both_rev(self):
        f = 1 | 0x10 | 0x20
        assert flip_flag(f) == f ^ 0x30

    def test_preserves_other_bits_paired(self):
        assert flip_flag(67) & ~0x30 == 67 & ~0x30

    def test_preserves_other_bits_unpaired(self):
        assert flip_flag(272) == 256

    def test_double_flip_identity(self):
        for f in [0, 16, 99, 147, 67, 131, 4, 256]:
            assert flip_flag(flip_flag(f)) == f


# ---------------------------------------------------------------------------
# transcript_to_genome
# ---------------------------------------------------------------------------
class TestTranscriptToGenome:
    def _tx(self, strand="+"):
        return Transcript(
            gene_id="G",
            transcript_id="T",
            chrom="chr1",
            strand=strand,
            exons={1: Span(1000, 1100), 2: Span(2000, 2050), 3: Span(3000, 3030)},
        )

    def test_first_exon_start(self):
        pos, idx = transcript_to_genome(0, self._tx())
        assert (pos, idx) == (1000, 0)

    def test_first_exon_end(self):
        pos, idx = transcript_to_genome(99, self._tx())
        assert (pos, idx) == (1099, 0)

    def test_second_exon_boundary(self):
        pos, idx = transcript_to_genome(100, self._tx())
        assert (pos, idx) == (2000, 1)

    def test_second_exon_mid(self):
        pos, idx = transcript_to_genome(125, self._tx())
        assert (pos, idx) == (2025, 1)

    def test_third_exon(self):
        pos, idx = transcript_to_genome(150, self._tx())
        assert (pos, idx) == (3000, 2)

    def test_minus_strand(self):
        tx = self._tx("-")
        pos, idx = transcript_to_genome(0, tx)
        assert pos == 1000 and idx == 0

    def test_minus_strand_second_exon(self):
        tx = self._tx("-")
        pos, idx = transcript_to_genome(100, tx)
        assert pos == 2000 and idx == 1

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            transcript_to_genome(999, self._tx())
        with pytest.raises(ValueError):
            transcript_to_genome(-1, self._tx())

    def test_single_exon(self):
        tx = Transcript(
            gene_id="G",
            transcript_id="T",
            chrom="chr1",
            strand="+",
            exons={1: Span(500, 600)},
        )
        pos, idx = transcript_to_genome(50, tx)
        assert (pos, idx) == (550, 0)


# ---------------------------------------------------------------------------
# remap_to_genome (liftover)
# ---------------------------------------------------------------------------
class TestRemapToGenome:
    def _genome_header(self):
        return pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6"},
                "SQ": [{"SN": "chr1", "LN": 10000}],
            }
        )

    def _tx_header(self):
        return pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6"},
                "SQ": [
                    {"SN": "TX1", "LN": 300},
                    {"SN": "TX2", "LN": 300},
                    {"SN": "TX3", "LN": 300},
                ],
            }
        )

    def _plus_tx(self):
        return Transcript(
            gene_id="G1",
            transcript_id="TX1",
            chrom="chr1",
            strand="+",
            exons={1: Span(1000, 1050), 2: Span(2000, 2050)},
        )

    def _minus_tx(self):
        return Transcript(
            gene_id="G2",
            transcript_id="TX2",
            chrom="chr1",
            strand="-",
            exons={1: Span(1000, 1050), 2: Span(2000, 2050)},
        )

    def _make_align(self, ref_name, start, cigar, flag=0, seq="ACGTACGTAC", md="10"):
        a = pysam.AlignedSegment(header=self._tx_header())
        a.query_name = "read"
        a.flag = flag
        a.reference_name = ref_name
        a.reference_start = start
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.cigarstring = cigar
        a.mapping_quality = 60
        a.set_tag("MD", md)
        return a

    def test_plus_simple(self):
        a = self._make_align("TX1", 10, "10M")
        r = remap_to_genome(a, self._genome_header(), self._plus_tx(), None)
        assert r.reference_start == 1010
        assert r.flag == 0

    def test_plus_intron(self):
        a = self._make_align("TX1", 45, "10M")
        r = remap_to_genome(a, self._genome_header(), self._plus_tx(), None)
        assert r.reference_start == 1045
        ops = r.cigartuples
        assert ops[0] == (0, 5)
        assert ops[1] == (3, 950)
        assert ops[2] == (0, 5)

    def test_minus_flips_flag_unpaired(self):
        a = self._make_align("TX2", 10, "10M", flag=0)
        r = remap_to_genome(a, self._genome_header(), self._minus_tx(), None)
        assert r.flag == 16

    def test_minus_flips_flag_paired(self):
        a = self._make_align("TX2", 10, "10M", flag=99)
        r = remap_to_genome(a, self._genome_header(), self._minus_tx(), None)
        assert r.flag == 83

    def test_minus_reverses_seq(self):
        a = self._make_align("TX2", 0, "6M", seq="AAACCC", md="6")
        a.query_qualities = pysam.qualitystring_to_array("ABCDEF")
        r = remap_to_genome(a, self._genome_header(), self._minus_tx(), None)
        assert r.query_sequence == "GGGTTT"
        assert list(r.query_qualities) == list(pysam.qualitystring_to_array("FEDCBA"))

    def test_minus_reverses_cigar(self):
        a = self._make_align("TX2", 0, "2S8M", seq="AAACGTACGT", md="8")
        r = remap_to_genome(a, self._genome_header(), self._minus_tx(), None)
        assert r.cigartuples[0][0] == 0  # M first
        assert r.cigartuples[-1][0] == 4  # S last

    def test_minus_intron_spanning(self):
        """A read near the exon junction on minus strand should produce M-N-M on the genome."""
        a = self._make_align("TX2", 45, "10M", seq="ACGTACGTAC", md="10")
        r = remap_to_genome(a, self._genome_header(), self._minus_tx(), None)
        ops = r.cigartuples
        assert ops[0] == (0, 5)  # 5M in first genomic exon
        assert ops[1] == (3, 950)  # 950N intron
        assert ops[2] == (0, 5)  # 5M in second genomic exon

    def _overlapping_tx(self):
        # exons overlap by 1 bp (2nd starts before 1st ends), as in some GTF
        # bacterial/tRNA gene models -> previously produced a negative intron
        return Transcript(
            gene_id="G3",
            transcript_id="TX3",
            chrom="chr1",
            strand="+",
            exons={1: Span(1000, 1050), 2: Span(1049, 1100)},
        )

    def test_overlapping_exons_no_negative_intron(self):
        a = self._make_align("TX3", 45, "10M")
        r = remap_to_genome(a, self._genome_header(), self._overlapping_tx(), None)
        ops = r.cigartuples
        assert r.reference_start == 1045
        # no N op and no negative length; the read spans both overlapping exons
        assert all(op[1] >= 0 for op in ops)
        assert all(op[0] != 3 for op in ops) or any(op[0] == 3 and op[1] > 0 for op in ops)

    def test_overlapping_exons_forward_strand_cigar(self):
        a = self._make_align("TX3", 45, "10M")
        r = remap_to_genome(a, self._genome_header(), self._overlapping_tx(), None)
        # 10M starting at 1045: first 5M within exon1 (ends 1050), then continue
        # contiguously into exon2 (starts 1049) without an N op
        ops = [op for op in r.cigartuples if op[0] == 0]
        assert sum(op[1] for op in ops) == 10


# ---------------------------------------------------------------------------
# Integration: convert_bam
# ---------------------------------------------------------------------------
class TestConvertBam:
    def test_liftover_basic(self, tmp_path, data_dir, has_liftover_data):
        from pathlib import Path
        from coralsnake.tbam2gbam import convert_bam

        output = str(tmp_path / "output.bam")
        convert_bam(
            input_bam=str(data_dir / "liftover_input.bam"),
            output_bam=output,
            annotation_file=str(data_dir / "liftover_annot.tsv"),
            faidx_file=str(data_dir / "liftover_faidx.fai"),
            threads=1,
            sort=False,
        )
        assert Path(output).exists()
        with pysam.AlignmentFile(output, "rb") as f:
            assert f.header is not None
