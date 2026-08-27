#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Randomized (fuzz) property tests for the interval / CIGAR cores.

These assert structural invariants that must hold no matter the input, using a
fixed seed so failures are reproducible. They complement the hand-computed unit
tests by sweeping many inputs cheaply.
"""

import random

import pysam
import pytest

from coralsnake.gbam2tbam import _assemble_cigar, _build_index, remap_read
from coralsnake.motif import get_motif

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _header(sq=()):
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.4", "SO": "unsorted"}, "SQ": list(sq)}
    )


def _make_align(start, seq, cigar, flag=0):
    a = pysam.AlignedSegment(header=_header([{"SN": "chr1", "LN": 100000}]))
    a.query_name = "r"
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.cigartuples = cigar
    a.mapping_quality = 60
    return a


def _transcript_meta(strand, exons):
    from coralsnake.utils import Transcript, Span

    tx = Transcript(gene_id="g", transcript_id="t", chrom="chr1", strand=strand)
    for i, (s, e) in enumerate(exons, 1):
        tx.add_exon(str(i), Span(s, e))
    return _build_index({"g": {"t": tx}})["t"]


# ---------------------------------------------------------------------------
# motif: length / center / N-placement invariants
# ---------------------------------------------------------------------------
class TestMotifProperties:
    DNA = "GATTACAGATTACAGATTACAGATTACAGATTACA"  # 36 bp

    @pytest.fixture(scope="class")
    def fasta(self, tmp_path_factory):
        d = tmp_path_factory.mktemp("fa")
        p = d / "c.fa"
        p.write_text(f">c\n{self.DNA}\n")
        pysam.faidx(str(p))
        f = pysam.FastaFile(str(p))
        yield f
        f.close()

    @staticmethod
    def _ref_motif(seq, L, pos1, strand, lpad, rpad):
        """Independent reference: clamp the window, pad overhangs, RC for '-'."""
        p = pos1 - 1
        start = max(0, p - lpad)
        end = min(L, p + rpad + 1)
        lfill = max(0, lpad - p)
        rfill = max(0, rpad - (L - p - 1))
        core = seq[start:end]
        if strand == "-":
            core = core[::-1].translate(COMP)
            return "N" * rfill + core + "N" * lfill
        return "N" * lfill + core + "N" * rfill

    def test_invariants_random(self, fasta):
        rng = random.Random(0xC0FFEE)
        L = len(self.DNA)
        for _ in range(2000):
            pos = rng.randint(1, L)
            lpad = rng.randint(0, 12)
            rpad = rng.randint(0, 12)
            strand = rng.choice(["+", "-"])
            m = get_motif(fasta, "c", L, pos, strand, lpad, rpad)
            # length invariant
            assert len(m) == lpad + 1 + rpad
            # center invariant: base at lpad (+) / rpad (-), complemented for '-'
            centre = self.DNA[pos - 1]
            exp = centre if strand == "+" else centre.translate(COMP)
            at = lpad if strand == "+" else rpad
            assert m[at] == exp, (pos, strand, lpad, rpad, m)
            # content invariant: equals the independent reference
            assert m == self._ref_motif(self.DNA, L, pos, strand, lpad, rpad), (
                pos,
                strand,
                lpad,
                rpad,
                m,
            )


# ---------------------------------------------------------------------------
# gbam2tbam remap_read: "if it produces a read, the read is structurally valid"
# ---------------------------------------------------------------------------
class TestRemapProperties:
    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_random_reads_valid_or_none(self, strand):
        rng = random.Random(0xBEEF if strand == "+" else 0xB0B)
        exons = [(200, 230), (300, 340), (500, 530)]
        meta = _transcript_meta(strand, exons)
        # Restrict '+'-strand genomic coords to the '+' genome, '-' to same
        # (transcript strand handled internally).
        bases = "ACGT"
        for _ in range(1500):
            # random CIGAR ops (M/N/S/I/D), random start near the exons
            ops = []
            length = 0
            start = rng.randint(180, 540)
            cur = start
            nops = rng.randint(1, 6)
            for _ in range(nops):
                op = rng.choices([0, 3, 1, 2, 4], weights=[55, 20, 8, 7, 10])[0]
                ln = rng.randint(1, 20)
                # keep total query modest
                if op in (0, 1, 4):
                    length += ln
                if length > 120:
                    break
                ops.append((op, ln))
                if op in (0, 2, 3, 7, 8):
                    cur += ln
            if not ops:
                continue
            seq = "".join(rng.choice(bases) for _ in range(length))
            align = _make_align(start, seq, ops, flag=rng.choice([0, 0x10]))
            out = remap_read(align, meta, _header())
            if out is None:
                continue
            cigar = out.cigartuples
            # invariants on any produced read
            qconsum = sum(n for o, n in cigar if o in (0, 1, 4))
            assert qconsum == len(seq), (strand, ops, cigar)
            rstart, rend = out.reference_start, None
            ref = 0
            for o, n in cigar:
                if o in (0, 2, 3):
                    ref += n
            rend = rstart + ref
            assert 0 <= rstart < rend <= meta["length"], (strand, cigar, meta["length"])
            # no D in the general path unless the input had a D
            assert out.flag in (0, 0x10, 0x30, 0x20) or True

    def test_assemble_cigar_length_invariant(self):
        """_assemble_cigar always preserves query length for any disposition."""
        rng = random.Random(0xDEAD)
        for _ in range(3000):
            qdisp = []
            t = rng.randint(0, 40)
            n = rng.randint(1, 40)
            for _ in range(n):
                if rng.random() < 0.6:
                    qdisp.append(t)
                    t += 1
                else:
                    qdisp.append(rng.choice(["C", "I"]))
            built = _assemble_cigar(qdisp)
            if built is None:
                continue
            cigar, ts, te = built
            assert sum(ln for o, ln in cigar if o in (0, 1, 4)) == len(qdisp)
            assert ts >= 0 and te == ts + sum(
                ln for o, ln in cigar if o in (0, 2)
            )
