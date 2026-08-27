#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for coralsnake.gbam2tbam (genome BAM -> transcript BAM remap)."""

from pathlib import Path

import pysam
import pytest

pytest.importorskip(
    "coralsnake.gbam2tbam", reason="gbam2tbam deps not available", exc_type=ImportError
)

from coralsnake.gbam2tbam import _build_index, _interval_to_transcript, remap_read  # noqa: E402


def _transcript(strand, spans):
    from coralsnake.utils import Transcript, Span

    tx = Transcript(gene_id="G", transcript_id="T", chrom="c", strand=strand)
    for i, (s, e) in enumerate(spans, 1):
        tx.add_exon(str(i), Span(s, e))
    return tx


def _header(sq=()):
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.4", "SO": "unsorted"}, "SQ": list(sq)}
    )


class TestBuildIndex:
    def test_plus_offsets(self):
        idx = _build_index({"g": {"t": _transcript("+", [(100, 110), (200, 230)])}})
        m = idx["t"]
        assert m["exons"] == [(100, 110), (200, 230)]
        assert m["offs"] == [0, 10]
        assert m["length"] == 40
        assert m["strand"] == "+"

    def test_minus_offsets_5p_first(self):
        # minus strand: 5' end is the rightmost exon
        idx = _build_index({"g": {"t": _transcript("-", [(100, 110), (200, 230)])}})
        m = idx["t"]
        assert m["exons"] == [(200, 230), (100, 110)]
        assert m["offs"] == [0, 30]
        assert m["length"] == 40


class TestIntervalToTranscript:
    def test_plus_inside_single_exon(self):
        idx = _build_index({"g": {"t": _transcript("+", [(100, 110), (200, 230)])}})
        m = idx["t"]
        assert _interval_to_transcript(102, 105, m) == (2, 5)   # exon1
        assert _interval_to_transcript(205, 208, m) == (15, 18)  # exon2 (offs+5)

    def test_minus_inside_single_exon(self):
        idx = _build_index({"g": {"t": _transcript("-", [(100, 110), (200, 230)])}})
        m = idx["t"]
        # 5' exon (genomic right, 200-230): tpos = 229 - g
        assert _interval_to_transcript(205, 208, m) == (22, 25)
        # 3' exon (genomic left, 100-110), offs=30: tpos = 30 + (109 - g)
        assert _interval_to_transcript(102, 105, m) == (35, 38)


class TestRemapRead:
    def _read(self, ref, start, seq, cigar, flag=0):
        a = pysam.AlignedSegment(header=_header([{"SN": "chr1", "LN": 1000}]))
        a.query_name = "r1"
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.cigartuples = cigar
        a.mapping_quality = 60
        return a

    def test_single_exon_read_plus(self):
        idx = _build_index({"g": {"t": _transcript("+", [(100, 110)])}})
        read = self._read("chr1", 102, "TTTTTTT", [(0, 7)])
        new = remap_read(read, idx["t"], _header())
        assert new.reference_start == 2
        assert new.cigartuples == [(0, 7)]
        assert new.flag == 0  # '+' transcript: flag unchanged

    def test_junction_read_joins_contiguous_on_transcript(self):
        idx = _build_index({"g": {"t": _transcript("+", [(100, 110), (200, 230)])}})
        # 8M in exon1 [102,110) -> t[2,10); intron to exon2 start; 3M [200,203) -> t[10,13)
        read = self._read("chr1", 102, "ACGTACGTACG", [(0, 8), (3, 90), (0, 3)])
        new = remap_read(read, idx["t"], _header())
        assert new.reference_start == 2
        assert new.cigartuples == [(0, 11)]  # contiguous on the transcript

    def test_skipped_exon_gets_ref_skip(self):
        # exons 100-110 and 200-230, read only in exon1 and exon3(200-230)
        idx = _build_index({"g": {"t": _transcript("+", [(100, 110), (200, 230)])}})
        # 5M at [102,107) -> t[2,7); 3M at [205,208) -> t[15,18) -> gap N(8)
        read = self._read("chr1", 102, "ACGTACGT", [(0, 5), (3, 98), (0, 3)])
        new = remap_read(read, idx["t"], _header())
        assert new.cigartuples == [(0, 5), (3, 8), (0, 3)]

    def test_minus_transcript_flips_strand_flag(self):
        idx = _build_index({"g": {"t": _transcript("-", [(100, 110), (200, 230)])}})
        read = self._read("chr1", 205, "ACGTACGT", [(0, 8)])  # forward on genome
        new = remap_read(read, idx["t"], _header())
        assert new.flag == 0x10  # forward genome read on '-' transcript => reverse

    def test_fully_intronic_read_returns_none(self):
        idx = _build_index({"g": {"t": _transcript("+", [(100, 110), (200, 230)])}})
        read = self._read("chr1", 150, "AAA", [(0, 3)])
        assert remap_read(read, idx["t"], _header()) is None

    def test_read_with_softclip_preserves_query_length(self):
        idx = _build_index({"g": {"t": _transcript("+", [(100, 120)])}})
        # 3S + 5M -> 8 query bases
        read = self._read("chr1", 105, "AAACCGGG", [(4, 3), (0, 5)])
        new = remap_read(read, idx["t"], _header())
        assert new.reference_start == 5
        assert new.cigartuples == [(4, 3), (0, 5)]

    def test_internal_deletion_read_skipped(self):
        idx = _build_index({"g": {"t": _transcript("+", [(100, 130)])}})
        read = self._read("chr1", 102, "ACGTACGT", [(0, 4), (2, 2), (0, 4)])
        assert remap_read(read, idx["t"], _header()) is None


class TestConvertEndToEnd:
    @pytest.fixture(params=[1, 4])
    def threads(self, request):
        return request.param

    def test_roundtrip_convert(self, tmp_path: Path, threads):
        from coralsnake.gbam2tbam import convert_bam

        annot = tmp_path / "annot.tsv"
        annot.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "g\tt\tchr1\t+\t100-109,200-229\n"
        )
        in_bam = tmp_path / "genome.bam"
        gh = _header([{"SN": "chr1", "LN": 1000}])
        with pysam.AlignmentFile(str(in_bam), "wb", header=gh) as out:
            a = pysam.AlignedSegment()
            a.query_name = "reads1"
            a.query_sequence = "ACGTACGTAC"  # 10 bases
            a.query_qualities = pysam.qualitystring_to_array("I" * 10)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 102
            # exon1 [99,109), intron, exon2 [199,229):
            # 7M [102,109) -> t[3,10); 90N -> g=199; 3M [199,202) -> t[10,13)
            a.cigartuples = [(0, 7), (3, 90), (0, 3)]
            a.mapping_quality = 60
            out.write(a)

        out_bam = tmp_path / "transcript.bam"
        convert_bam(str(in_bam), str(out_bam), str(annot), threads=threads)
        with pysam.AlignmentFile(str(out_bam), "rb") as r:
            reads = list(r)
        assert len(reads) == 1
        read = reads[0]
        assert read.reference_name == "t"
        assert read.reference_start == 3
        assert read.cigartuples == [(0, 10)]  # contiguous on the transcript
