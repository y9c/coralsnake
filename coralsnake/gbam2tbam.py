#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# gbam2tbam: remap a GENOME-aligned BAM back onto TRANSCRIPT references.
#
# Inverse of `tbam2gbam` (`coralsnake liftover`, `coralsnake tbam2gbam`).
# Reads aligned to genomic coordinates are clipped at exon boundaries and
# re-mapped to per-transcript reference sequences in 5'->3' orientation
# (matching the transcript FASTA emitted by `prepare`). Introns disappear, so a
# spliced read spanning an intron becomes contiguous on the transcript; a read
# that skips an exon instead gets a ref-skip (N) of the skipped exonic bases.
#
# Converts reads whose CIGAR is M / N / leading-trailing S (the canonical
# spliced transcriptome case). Reads with internal insertions or deletions, or
# with M blocks that extend into introns, are silently skipped rather than
# emitting an invalid/truncated alignment (documented limitation of v1).

import os
from functools import lru_cache

import pysam
from rich.progress import track

from .utils import get_logger, load_annotation

LOGGER = get_logger(__name__)


@lru_cache(maxsize=10000)
def flip_flag(flag):
    """Flip the read strand flag (0x10/0x20) for a '-' transcript."""
    if flag & 1:  # paired: toggle read-reverse (16) + mate-reverse (32)
        return flag ^ 0x30
    return flag ^ 16  # unpaired: toggle read-reverse (16)


def _exons_53(transcript):
    """Exons ordered 5'->3' plus each exon's 5' transcript offset and total length."""
    exons = sorted(transcript.exons.values(), key=lambda e: e.start)
    if transcript.strand == "-":
        exons.reverse()  # for '-', the 5' end is the rightmost exon
    spans = [(e.start, e.end) for e in exons]
    offs = []
    acc = 0
    for gs, ge in spans:
        offs.append(acc)
        acc += ge - gs
    return spans, offs, acc


def _build_index(annot):
    """Flatten {gene_id: {tid: Transcript}} into {tid: index_dict}."""
    flat = {}

    def add(tid, tx):
        spans, offs, tlen = _exons_53(tx)
        flat[tid] = {
            "exons": spans,
            "offs": offs,
            "length": tlen,
            "strand": tx.strand,
        }

    for g_id in annot:
        for t_id, tx in annot[g_id].items():
            add(t_id, tx)
        if len(annot[g_id]) == 1:
            add(g_id, next(iter(annot[g_id].values())))
    return flat


def _interval_to_transcript(gs, ge, meta):
    """Map a fully-exonic genomic interval [gs,ge) -> transcript [ts,te). None if not exonic."""
    ts = te = None
    for i, (es, ee) in enumerate(meta["exons"]):
        lo = max(gs, es)
        hi = min(ge, ee)
        if lo >= hi:
            continue
        if meta["strand"] == "+":
            a = meta["offs"][i] + (lo - es)
            b = meta["offs"][i] + (hi - es)
        else:
            a = meta["offs"][i] + (ee - hi)
            b = meta["offs"][i] + (ee - lo)
        ts = a if ts is None else min(ts, a)
        te = b if te is None else max(te, b)
    if ts is None:
        return None
    return ts, te


def _fully_exonic(meta, gs, ge):
    """True if every base of [gs,ge) lies inside the transcript's exons."""
    covered = 0
    for es, ee in meta["exons"]:
        covered += max(0, min(ge, ee) - max(gs, es))
    return covered == ge - gs


def _classify_cigar(cigartuples):
    """Return (leading_s, trailing_s, middle) where middle ops must be only M/N.

    Returns (leading_s, trailing_s, middle_ops, clean). ``middle_ops`` is the
    list of (op, len) ops strictly between the leading/trailing soft-clips.
    ``clean`` is False if the middle contains anything other than M(0)/N(3).
    """
    tuples = list(cigartuples or [])
    leading = trailing = 0
    # leading soft-clips
    i = 0
    while i < len(tuples) and tuples[i][0] == 4:
        leading += tuples[i][1]
        i += 1
    # trailing soft-clips
    j = len(tuples)
    while j > i and tuples[j - 1][0] == 4:
        trailing += tuples[j - 1][1]
        j -= 1
    middle = tuples[i:j]
    clean = all(op in (0, 3) for op, _ in middle)
    return leading, trailing, middle, clean


def remap_read(align, meta, t_header):
    """Remap one genome-aligned read onto a transcript reference (5'->3').

    Returns a new AlignedSegment, or None if the read cannot be cleanly mapped.
    """
    if align.cigartuples is None:
        return None
    leading_s, trailing_s, middle, clean = _classify_cigar(align.cigartuples)
    if not clean:
        return None

    # Aligned reference blocks (M/=/X), excluding N introns. Verify each is
    # fully exonic; map to transcript coordinates.
    mapped = []
    for (gs, ge) in align.get_blocks():
        if not _fully_exonic(meta, gs, ge):
            return None
        iv = _interval_to_transcript(gs, ge, meta)
        if iv is None:
            return None
        mapped.append(iv)
    if not mapped:
        return None
    mapped.sort()
    merged = [list(mapped[0])]
    for a, b in mapped[1:]:
        if a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])

    t_start = merged[0][0]
    cigar = []
    prev_end = None
    for a, b in merged:
        if prev_end is not None and a > prev_end:
            cigar.append((3, a - prev_end))  # skipped exonic gap on the transcript
        cigar.append((0, b - a))
        prev_end = b
    if leading_s:
        cigar.insert(0, (4, leading_s))
    if trailing_s:
        cigar.append((4, trailing_s))

    # Safety: the query- (M+I+S) consuming length must equal the read length.
    if sum(length for op, length in cigar if op in (0, 1, 4)) != len(align.query_sequence):
        return None

    new = pysam.AlignedSegment(header=t_header)
    new.query_name = align.query_name
    new.query_sequence = align.query_sequence
    new.query_qualities = align.query_qualities
    # A '-' transcript is the reverse complement of the genomic '+' strand, so a
    # read forward on the genome is reverse on the transcript (and vice-versa).
    new.flag = flip_flag(align.flag) if meta["strand"] == "-" else align.flag
    for key, value in align.get_tags():
        if key in ("NM", "MD", "AS", "XS", "YS", "XA"):  # reference-position-dependent
            continue
        new.set_tag(key, value)
    new.reference_start = t_start
    new.cigartuples = cigar
    new.mapping_quality = align.mapping_quality
    return new


def _overlaps(meta, align):
    rs, re = align.reference_start, align.reference_end
    return any(es < re and ee > rs for (es, ee) in meta["exons"])


def _overlap_len(meta, align):
    rs, re = align.reference_start, align.reference_end
    return sum(max(0, min(re, ee) - max(rs, es)) for (es, ee) in meta["exons"])


def _pick_transcript(annot, align):
    """Choose the best transcript: prefer full containment, then most overlap."""
    best, best_score = None, -1
    for tid, meta in annot.items():
        if not _overlaps(meta, align):
            continue
        contained = all(
            any(es <= b0 and ee >= b1 for (es, ee) in meta["exons"])
            for (b0, b1) in align.get_blocks()
        )
        score = (10**9 if contained else 0) + _overlap_len(meta, align)
        if score > best_score:
            best, best_score = tid, score
    return best


def _sort_out_bam(output_bam, threads):
    if not output_bam.endswith(".bam"):
        return
    LOGGER.info("Sorting output BAM file...")
    import shutil
    import tempfile

    out_dir = os.path.dirname(os.path.abspath(output_bam))
    fd, tmp = tempfile.mkstemp(dir=out_dir, prefix=".coralsnake_sort_", suffix=".bam")
    os.close(fd)
    os.remove(tmp)
    try:
        pysam.sort("-@", str(threads), "-o", tmp, output_bam)
        shutil.move(tmp, output_bam)
        pysam.index("-@", str(threads), output_bam)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def convert_bam(input_bam, output_bam, annotation_file, threads=8, sort=False):
    """Remap a genome-aligned BAM to transcript references (gbam -> tbam)."""
    LOGGER.info("Loading annotation...")
    annot = _build_index(load_annotation(annotation_file))

    t_dict = {
        "HD": {"VN": "1.4", "SO": "unsorted"},
        "SQ": [{"SN": tid, "LN": meta["length"]} for tid, meta in annot.items()],
    }
    t_header = pysam.AlignmentHeader.from_dict(t_dict)

    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        with pysam.AlignmentFile(
            output_bam, "wb" if output_bam.endswith(".bam") else "w", header=t_header
        ) as out_bam:
            for align in track(in_bam, description="Processing..."):
                if align.is_unmapped or align.reference_name is None:
                    continue
                tid = _pick_transcript(annot, align)
                if tid is None:
                    continue
                new = remap_read(align, annot[tid], t_header)
                if new is None:
                    continue
                new.reference_name = tid
                out_bam.write(new)

    if sort:
        _sort_out_bam(output_bam, threads)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Remap genome BAM to transcript BAM")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-a", "--annotation", required=True)
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("-s", "--sort", action="store_true")
    args = parser.parse_args()
    convert_bam(args.input, args.output, args.annotation, args.threads, args.sort)
