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


def _per_query_disposition(align, meta):
    """Walk a genome CIGAR into a per-query-base transcript disposition.

    Returns ``(qdisp, qlen, has_deletion)``, where ``qdisp`` is a list parallel
    to the read's query bases whose entries are an int (transcript position for
    an exonic M/=/X base), ``'C'`` (soft-clip: explicit soft-clip, or an M base
    that fell inside an intron), or ``'I'`` (insertion). Returns ``None`` if the
    read cannot be cleanly handled (e.g. any reference deletion ``D``, which is
    not yet remapped).
    """
    if align.cigartuples is None:
        return None
    cigar = align.cigartuples
    if any(op == 2 for op, _ in cigar):
        return None  # deletions unsupported (conservative, rare)
    strand = meta["strand"]
    exons = meta["exons"]
    offs = meta["offs"]

    def tpos(idx, g):
        es, ee = exons[idx]
        return offs[idx] + (g - es if strand == "+" else ee - 1 - g)

    def containing(g):
        for i, (es, ee) in enumerate(exons):
            if es <= g < ee:
                return i
        return None

    seq = align.query_sequence or ""
    qlen = len(seq)
    qdisp = []
    g = align.reference_start
    q = 0
    for op, length in cigar:
        if op in (0, 7, 8):  # M/=/X: query + reference
            for _ in range(length):
                idx = containing(g)
                qdisp.append(tpos(idx, g) if idx is not None else "C")
                q += 1
                g += 1
        elif op == 1:  # I: query only
            for _ in range(length):
                qdisp.append("I")
                q += 1
        elif op == 3:  # N: reference-only (intron) - advance genome
            g += length
        elif op == 4:  # S: query only (soft clip)
            for _ in range(length):
                qdisp.append("C")
                q += 1
        elif op in (5, 8):  # H, P
            pass
    if q != qlen:
        return None  # input CIGAR/query length mismatch
    return qdisp, qlen, False


def _assemble_cigar(qdisp):
    """Build a transcript CIGAR from per-query transcript dispositions.

    ``qdisp`` must already be in transcript (5'->3', increasing position) order.
    Returns ``(cigar_tuples, reference_start, reference_end)`` or ``None``.
    """
    out = []
    last_t = None
    first_t = None
    for it in qdisp:
        if isinstance(it, int):
            t = it
            if first_t is None:
                first_t = t
            if last_t is None:
                out.append([0, 1])
            elif t == last_t + 1:
                if out and out[-1][0] == 0:
                    out[-1][1] += 1
                else:
                    out.append([0, 1])
            elif t > last_t + 1:
                out.append([3, t - last_t - 1])  # skipped exonic gap
                out.append([0, 1])
            else:
                return None  # overlapping bases - bail
            last_t = t
        elif it == "I":
            if out and out[-1][0] == 1:
                out[-1][1] += 1
            else:
                out.append([1, 1])
        elif it == "C":
            if out and out[-1][0] == 4:
                out[-1][1] += 1
            else:
                out.append([4, 1])
    if first_t is None:
        return None
    return [(op, n) for op, n in out], first_t, last_t + 1


def remap_read(align, meta, t_header):
    """Remap one genome-aligned read onto a transcript reference (5'->3').

    Handles M/=/X, N (intron), soft-clips, insertions, and M blocks that dip
    into introns (those bases become soft-clips on the transcript). Returns a
    new AlignedSegment, or None if the read cannot be cleanly mapped.
    """
    res = _per_query_disposition(align, meta)
    if res is None:
        return None
    qdisp, _qlen, _has_del = res
    strand = meta["strand"]

    # For '+' transcripts the read's query order equals transcript (increasing
    # position) order. For '-' transcripts it is reversed (the read is written
    # reverse, flag flipped), so iterate the disposition in reverse.
    qseq = qdisp if strand == "+" else qdisp[::-1]
    built = _assemble_cigar(qseq)
    if built is None:
        return None
    cigar, t_start, t_end = built

    seq = align.query_sequence or ""
    if sum(n for op, n in cigar if op in (0, 1, 4)) != len(seq):
        return None
    if t_start < 0 or t_end > meta["length"]:
        return None

    new = pysam.AlignedSegment(header=t_header)
    new.query_name = align.query_name
    new.query_sequence = align.query_sequence
    new.query_qualities = align.query_qualities
    # A '-' transcript is the reverse complement of the '+' genomic forward
    # strand, so a read forward on the genome is reverse on the transcript
    # (and vice-versa): flip the strand flag.
    new.flag = flip_flag(align.flag) if strand == "-" else align.flag
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
            if threads <= 1:
                _convert_single(in_bam, out_bam, annot, t_header)
            else:
                _convert_parallel(in_bam, out_bam, annot, t_header, threads)

    if sort:
        _sort_out_bam(output_bam, threads)


def _remap_one(align, annot, t_header):
    """Remap a single read; return an AlignedSegment or None (skipped)."""
    if align.is_unmapped or align.reference_name is None:
        return None
    tid = _pick_transcript(annot, align)
    if tid is None:
        return None
    new = remap_read(align, annot[tid], t_header)
    if new is None:
        return None
    new.reference_name = tid
    return new


def _convert_single(in_bam, out_bam, annot, t_header):
    for align in track(in_bam, description="Processing..."):
        new = _remap_one(align, annot, t_header)
        if new is not None:
            out_bam.write(new)


# Worker globals for multiprocessing
_WORKER_ANNOT = None
_WORKER_THEADER = None
_WORKER_GHEADER = None


def _init_worker(annot_flattened, genome_header_dict, t_header_dict):
    global _WORKER_ANNOT, _WORKER_THEADER, _WORKER_GHEADER
    _WORKER_ANNOT = annot_flattened
    _WORKER_THEADER = pysam.AlignmentHeader.from_dict(t_header_dict)
    _WORKER_GHEADER = pysam.AlignmentHeader.from_dict(genome_header_dict)


def _process_chunk(chunk_strings):
    """Remap a chunk of serialised genome reads; return serialised transcript reads."""
    out = []
    # Input reads are GENOME-aligned -> deserialize with the genome header so
    # their reference_name resolves (the transcript header only knows transcripts).
    gh = _WORKER_GHEADER
    th = _WORKER_THEADER
    for s in chunk_strings:
        align = pysam.AlignedSegment.fromstring(s, gh)
        new = _remap_one(align, _WORKER_ANNOT, th)
        if new is not None:
            out.append(new.to_string())
    return out


def _convert_parallel(in_bam, out_bam, annot, t_header, threads):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    t_dict = t_header.to_dict()
    genome_dict = in_bam.header.to_dict()
    chunk_size = 5000
    max_queue = threads * 2
    with ProcessPoolExecutor(
        max_workers=threads,
        mp_context=mp.get_context("spawn"),
        initializer=_init_worker,
        initargs=(annot, genome_dict, t_dict),
    ) as executor:
        futures = []
        chunk = []
        for align in track(in_bam, description="Processing..."):
            chunk.append(align.to_string())
            if len(chunk) >= chunk_size:
                futures.append(executor.submit(_process_chunk, list(chunk)))
                chunk = []
                if len(futures) >= max_queue:
                    _drain_one(futures, t_header, out_bam)
        if chunk:
            futures.append(executor.submit(_process_chunk, chunk))
        while futures:
            _drain_one(futures, t_header, out_bam)


def _drain_one(futures, t_header, out_bam):
    for s in futures.pop(0).result():
        out_bam.write(pysam.AlignedSegment.fromstring(s, t_header))


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
