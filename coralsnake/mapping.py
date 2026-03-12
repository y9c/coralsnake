#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-08 20:32


import atexit
import multiprocessing as mp
import os
import random
import re
import shutil
import sys
import tempfile
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from contextlib import ExitStack
from functools import lru_cache

import pysam
from bwamem import BwaAligner, BwaIndexer, fastx_read, read_paired_fastx
from bwamem.libbwa import ffi, libbwa
from rich.progress import Progress, SpinnerColumn, TextColumn

from . import seqops
from .utils import (
    convert_file_realtime,
    format_duration,
    km_conversion,
    mk_conversion,
    reverse_complement,
)


@lru_cache(maxsize=100)
def _get_ref_cached(idx_ptr, rid, start, end):
    res_ptr = libbwa.bwa_fetch_seq(idx_ptr, rid, start, end)
    if res_ptr == ffi.NULL:
        return None
    res = ffi.string(res_ptr).decode()
    libbwa.free(res_ptr)
    return res


def _build_indices_with_progress(
    ref_file: str,
    index_base_dir: str,
    on_update,
    poll_interval: float = 0.1,
    ref_suffix: str = "",
):
    from concurrent.futures import ThreadPoolExecutor

    os.makedirs(index_base_dir, exist_ok=True)
    orig_fa = os.path.join(index_base_dir, f"ref{ref_suffix}.orig.fa")
    orig_prefix = os.path.splitext(orig_fa)[0]
    mk_fa = os.path.join(index_base_dir, f"ref{ref_suffix}.mk.fa")
    mk_prefix = os.path.splitext(mk_fa)[0]

    if ref_file is None:
        on_update("Done", 100.0, 100.0)
        return orig_fa, mk_prefix

    if os.path.abspath(ref_file) != os.path.abspath(orig_fa):
        shutil.copyfile(ref_file, orig_fa)

    indexer1 = BwaIndexer(verbose=2, capture_progress=True)
    indexer2 = BwaIndexer(verbose=2, capture_progress=True)

    def build_orig():
        indexer1.build_index(orig_fa, prefix=orig_prefix)

    def build_mk():
        indexer2.build_index(mk_fa, prefix=mk_prefix)

    def do_conversion():
        convert_file_realtime(ref_file, mk_fa, "AC", "GT")

    with ThreadPoolExecutor(max_workers=3) as ex:
        fut_orig = ex.submit(build_orig)
        fut_conv = ex.submit(do_conversion)
        while not fut_conv.done():
            p1 = 100.0 if fut_orig.done() else (indexer1.progress_percent or 0.0)
            on_update("Converting...", p1, 0.0)
            time.sleep(poll_interval)
        fut_conv.result()
        fut_mk = ex.submit(build_mk)
        while True:
            p1 = 100.0 if fut_orig.done() else (indexer1.progress_percent or 0.0)
            p2 = 100.0 if fut_mk.done() else (indexer2.progress_percent or 0.0)
            on_update("Indexing...", p1, p2)
            if fut_orig.done() and fut_mk.done():
                break
            time.sleep(poll_interval)
        fut_orig.result()
        fut_mk.result()

    on_update("Done", 100.0, 100.0)
    return orig_fa, mk_prefix


def _map_batch_worker(
    batch, paired, max_mismatches, min_alignment_length, min_mapping_ratio
):
    results = []
    if not paired:
        seqs = [item[1] for item in batch]
        conv1 = (
            seqops.batch_base_conversion(seqs, "AC", "GT")
            if _FORWARD_LIBRARY
            else seqops.batch_base_conversion(seqs, "GT", "AC")
        )
        conv2 = (
            seqops.batch_base_conversion(seqs, "GT", "AC")
            if _FORWARD_LIBRARY
            else seqops.batch_base_conversion(seqs, "AC", "GT")
        )
        for i, (name1, seq1, qua1) in enumerate(batch):
            ref_idx, mapping_result = run_mapping_se(
                seq1,
                qua1,
                conv1[i],
                conv2[i],
                max_mismatches,
                min_alignment_length,
                min_mapping_ratio,
            )
            results.append((ref_idx, mapping_result))
    else:
        seqs1 = [item[0][1] for item in batch]
        seqs2 = [item[1][1] for item in batch]
        rc_s1 = [seqops.reverse_complement(s) for s in seqs1]
        rc_s2 = [seqops.reverse_complement(s) for s in seqs2]
        if _FORWARD_LIBRARY:
            c1_r1, c1_r2 = (
                seqops.batch_base_conversion(seqs1, "AC", "GT"),
                seqops.batch_base_conversion(rc_s2, "AC", "GT"),
            )
            c2_r1, c2_r2 = (
                seqops.batch_base_conversion(rc_s1, "AC", "GT"),
                seqops.batch_base_conversion(seqs2, "AC", "GT"),
            )
        else:
            c1_r1, c1_r2 = (
                seqops.batch_base_conversion(rc_s1, "AC", "GT"),
                seqops.batch_base_conversion(seqs2, "AC", "GT"),
            )
            c2_r1, c2_r2 = (
                seqops.batch_base_conversion(seqs1, "AC", "GT"),
                seqops.batch_base_conversion(rc_s2, "AC", "GT"),
            )
        conv_data = (c1_r1, c1_r2, c2_r1, c2_r2)
        for i, ((name1, seq1, qua1), (name2, seq2, qua2)) in enumerate(batch):
            ref_idx, mapping_result = run_mapping_pe(
                seq1,
                seq2,
                qua1,
                qua2,
                i,
                conv_data,
                rc_s1[i],
                rc_s2[i],
                max_mismatches,
                min_alignment_length,
                min_mapping_ratio,
            )
            results.append((ref_idx, mapping_result))
    return results


def _check_reference_length(aligner, min_length=100):
    for i in range(aligner.index.bns.n_seqs):
        if aligner.index.bns.anns[i].len >= min_length:
            return True
    return False


def run_mapping_se(
    seq1, qua1, s1_c1, s1_c2, max_mismatches, min_alignment_length, min_mapping_ratio
):
    rc_seq1 = None
    for ref_idx in range(len(_ALIGNERS_MK)):
        if not _check_reference_length(_ALIGNERS_MK[ref_idx]):
            continue
        mapped = []
        oris = [1, 2] if _ORIENTATION_FILTER is None else [_ORIENTATION_FILTER]
        for orientation in oris:
            s1_c = s1_c1 if orientation == 1 else s1_c2
            hits = _ALIGNERS_MK[ref_idx].align(s1_c, min_mapq=1)
            if not hits:
                continue
            for h in hits:
                if (h[2] - h[1]) < min_alignment_length:
                    continue
                is_rev = h[3] == -1
                s = seqops.reverse_complement(seq1) if is_rev else seq1
                rid = h[10] if len(h) > 10 else _RID_MAPS[ref_idx].get(h[0], -1)
                if rid < 0:
                    continue
                ref = _get_ref_cached(
                    _ALIGNERS_ORIG[ref_idx].index, rid, h[1], h[2] + 100
                )
                res = seqops.score_and_tag(h[7], s, ref, (orientation == 1) ^ is_rev)
                if (
                    not res
                    or res[0] < -500
                    or res[1] > max_mismatches
                    or (h[5] - h[4]) < len(seq1) * min_mapping_ratio
                ):
                    continue
                # SE data: [is_rev, rname, pos, mq, cigar, score, md, ori, yf, zf, yc, zc, ns, nc, rid]
                mapped.append(
                    [
                        is_rev,
                        h[0],
                        h[1],
                        score_to_mapq(res[0]),
                        h[7],
                        res[0],
                        res[2],
                        orientation,
                        res[3],
                        res[4],
                        res[5],
                        res[6],
                        res[7],
                        res[8],
                        rid,
                    ]
                )
        if mapped:
            mapped.sort(key=lambda x: x[5], reverse=True)
            return (ref_idx, mapped)
    return (None, [])


def run_mapping_pe(
    seq1,
    seq2,
    qua1,
    qua2,
    idx,
    conv_data,
    rc_s1,
    rc_s2,
    max_mismatches,
    min_alignment_length,
    min_mapping_ratio,
):
    q1_l, q2_l = len(seq1), len(seq2)
    c1_r1, c1_r2, c2_r1, c2_r2 = conv_data
    all_res = []
    for ref_idx in range(len(_ALIGNERS_MK)):
        if not _check_reference_length(_ALIGNERS_MK[ref_idx]):
            continue
        mapped = []
        oris = [1, 2] if _ORIENTATION_FILTER is None else [_ORIENTATION_FILTER]
        for orientation in oris:
            is_o1 = orientation == 1
            s1_c, s2_c = (c1_r1[idx], c1_r2[idx]) if is_o1 else (c2_r1[idx], c2_r2[idx])
            hits = _ALIGNERS_MK[ref_idx].align(s1_c, s2_c, min_mapq=1)
            if not hits:
                continue
            for h1, h2, is_p, isize in hits:
                res1 = res2 = None
                rid1 = rid2 = -1

                # Strand logic strictly matching standard FR orientation
                # If Library is Forward, O1 is forward, Read 1 maps Fwd, Read 2 maps Rev
                # The BAM flag `is_reverse` should track the biological strand of the fragment.
                if is_o1:
                    r1_abs_rev = (h1[3] == -1) if h1 else False
                    # Read 2 in O1 should map reverse if proper, so its true strand is reversed from its mapped strand
                    r2_abs_rev = (h2[3] == 1) if h2 else False
                else:
                    # In O2, the original fragment was reverse.
                    r1_abs_rev = (h1[3] == 1) if h1 else False
                    r2_abs_rev = (h2[3] == -1) if h2 else False

                if h1 and (h1[2] - h1[1]) >= min_alignment_length:
                    r1_rev = (
                        (h1[3] == -1 if _FORWARD_LIBRARY else h1[3] == 1)
                        if is_o1
                        else (h1[3] == 1 if _FORWARD_LIBRARY else h1[3] == -1)
                    )
                    s1_cur = rc_s1 if r1_rev else seq1
                    rid1 = h1[10] if len(h1) > 10 else _RID_MAPS[ref_idx].get(h1[0], -1)
                    if rid1 >= 0:
                        ref1 = _get_ref_cached(
                            _ALIGNERS_ORIG[ref_idx].index, rid1, h1[1], h1[2] + 100
                        )
                        res1 = seqops.score_and_tag(
                            h1[7],
                            s1_cur,
                            ref1,
                            (is_o1 if _FORWARD_LIBRARY else (not is_o1)) ^ r1_rev,
                        )
                if h2 and (h2[2] - h2[1]) >= min_alignment_length:
                    r2_rev = (
                        (h2[3] == 1 if _FORWARD_LIBRARY else h2[3] == -1)
                        if is_o1
                        else (h2[3] == -1 if _FORWARD_LIBRARY else h2[3] == 1)
                    )
                    s2_cur = rc_s2 if r2_rev else seq2
                    rid2 = h2[10] if len(h2) > 10 else _RID_MAPS[ref_idx].get(h2[0], -1)
                    if rid2 >= 0:
                        ref2 = _get_ref_cached(
                            _ALIGNERS_ORIG[ref_idx].index, rid2, h2[1], h2[2] + 100
                        )
                        res2 = seqops.score_and_tag(
                            h2[7],
                            s2_cur,
                            ref2,
                            (not (is_o1 if _FORWARD_LIBRARY else (not is_o1))) ^ r2_rev,
                        )

                if res1 or res2:
                    if (res1 and res1[0] < -500) or (res2 and res2[0] < -500):
                        continue
                    if (
                        (res1[1] if res1 else 0) + (res2[1] if res2 else 0)
                    ) > max_mismatches:
                        continue
                    if h1 and (h1[5] - h1[4]) < q1_l * min_mapping_ratio:
                        continue
                    if h2 and (h2[5] - h2[4]) < q2_l * min_mapping_ratio:
                        continue
                    mq = score_to_mapq(
                        min(res1[0] if res1 else 99, res2[0] if res2 else 99)
                    )
                    # PE data: [is_rev, rname, pos, mq, cigar, rnext, pnext, isize, score, md, ori, yf, zf, yc, zc, ns, nc, rid, is_p, mate_rev, mate_unmap]
                    m1 = (
                        [
                            r1_abs_rev,
                            h1[0],
                            h1[1],
                            mq,
                            h1[7],
                            h2[0] if h2 else "*",
                            h2[1] if h2 else 0,
                            (isize if h1[1] <= h2[1] else -isize) if h2 else 0,
                            res1[0],
                            res1[2],
                            orientation,
                            res1[3],
                            res1[4],
                            res1[5],
                            res1[6],
                            res1[7],
                            res1[8],
                            rid1,
                            bool(is_p),
                            bool(r2_abs_rev),
                            bool(not h2),
                        ]
                        if res1
                        else None
                    )
                    m2 = (
                        [
                            r2_abs_rev,
                            h2[0],
                            h2[1],
                            mq,
                            h2[7],
                            h1[0] if h1 else "*",
                            h1[1] if h1 else 0,
                            (isize if h2[1] <= h1[1] else -isize) if h1 else 0,
                            res2[0],
                            res2[2],
                            orientation,
                            res2[3],
                            res2[4],
                            res2[5],
                            res2[6],
                            res2[7],
                            res2[8],
                            rid2,
                            bool(is_p),
                            bool(r1_abs_rev),
                            bool(not h1),
                        ]
                        if res2
                        else None
                    )
                    mapped.append(
                        [(res1[0] if res1 else 0) + (res2[0] if res2 else 0), m1, m2]
                    )
        if mapped:
            mapped.sort(key=lambda x: x[0], reverse=True)
            if any(r[1] and r[2] for r in mapped):
                return (ref_idx, mapped)
            all_res.append((ref_idx, mapped))
        else:
            all_res.append((ref_idx, []))
    for r_idx, rs in all_res:
        if rs:
            return (r_idx, rs)
    return (None, [])


_ALIGNERS_ORIG = []
_ALIGNERS_MK = []
_ORIENTATION_FILTER = None
_FORWARD_LIBRARY = None
_RID_MAPS = []


def _init_worker(ref_indices, orientation_filter, forward_library):
    global \
        _ALIGNERS_ORIG, \
        _ALIGNERS_MK, \
        _ORIENTATION_FILTER, \
        _FORWARD_LIBRARY, \
        _RID_MAPS
    _ORIENTATION_FILTER, _FORWARD_LIBRARY = orientation_filter, forward_library
    _ALIGNERS_ORIG = []
    _ALIGNERS_MK = []
    _RID_MAPS = []
    for o_fa, mk_pre, _ in ref_indices:
        opts = {
            "min_seed_len": 14,
            "max_occ": 1000,
            "softclip_supplementary": True,
            "mark_secondary": True,
            "clip_penalties": (6, 6),
            "unpaired_penalty": 24,
            "min_score": 20,
            "insert_model": (80, 60, 450),
        }
        a_orig = BwaAligner(os.path.splitext(o_fa)[0], **opts)
        a_mk = BwaAligner(mk_pre, **opts)
        _ALIGNERS_ORIG.append(a_orig)
        _ALIGNERS_MK.append(a_mk)
        _RID_MAPS.append(
            {
                ffi.string(a_orig.index.bns.anns[i].name).decode(): i
                for i in range(a_orig.index.bns.n_seqs)
            }
        )
        import ctypes

        ptr = int(ffi.cast("uintptr_t", a_orig.index))


def _build_and_check_indices(ref_files, ref_indices, index_dirs, index_only):
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold green]Index[/bold green]"),
        TextColumn("{task.description}"),
        TextColumn(
            "| [yellow]{task.fields[status]}[/yellow] | [cyan]{task.fields[progress]}[/cyan]"
        ),
        transient=False,
    ) as progress:
        task = progress.add_task(
            f"Indices...", total=None, status="Checking", progress="0%"
        )
        for i in range(len(ref_indices)):
            ref_suffix = ref_indices[i][2]
            d = index_dirs[i] if len(index_dirs) == len(ref_indices) else index_dirs[0]
            o_fa, mk_pre = (
                os.path.join(d, f"ref{ref_suffix}.orig.fa"),
                os.path.join(d, f"ref{ref_suffix}.mk"),
            )
            if index_only or not (
                os.path.exists(o_fa)
                and all(
                    os.path.exists(mk_pre + e)
                    for e in [".amb", ".ann", ".bwt", ".pac", ".sa"]
                )
            ):

                def on_up(s, p1, p2):
                    progress.update(
                        task,
                        status=f"Ref {i + 1}: {s}",
                        progress=f"{(p1 + p2) / 2:.1f}%",
                    )

                idx0, idx_mk = _build_indices_with_progress(
                    ref_files[i] if i < len(ref_files) else None,
                    d,
                    on_up,
                    ref_suffix=ref_suffix,
                )
                ref_indices[i] = (idx0, idx_mk, ref_suffix)
            else:
                ref_indices[i] = (o_fa, mk_pre, ref_suffix)
        progress.update(
            task, description="✓ Indices ready", status="Done", progress="100%"
        )
        return ref_indices


def score_to_mapq(score):
    return max(0, min(60, score))


def create_bam_record(
    header, map_data, name, seq, qual, is_secondary, global_rid_map, read_id
):
    a = pysam.AlignedSegment(header=header)
    a.query_name = name
    a.is_secondary = is_secondary
    a.reference_id = global_rid_map.get(map_data[1], -1)
    a.reference_start, a.mapping_quality, a.cigarstring = (
        map_data[2],
        map_data[3],
        map_data[4],
    )
    a.is_reverse = bool(map_data[0])

    if len(map_data) > 15:  # PE Hit
        a.is_paired = True
        a.is_read1 = read_id == 1
        a.is_read2 = read_id == 2
        a.is_proper_pair = bool(map_data[18])
        a.mate_is_reverse = bool(map_data[19])
        a.mate_is_unmapped = bool(map_data[20])
        a.next_reference_id = global_rid_map.get(map_data[5], -1)
        a.next_reference_start = map_data[6] if not a.mate_is_unmapped else 0
        a.template_length = map_data[7]
        off = 8
    else:  # SE Hit
        a.is_paired = False
        off = 5

    if a.is_reverse:
        a.query_sequence = reverse_complement(seq)
        a.query_qualities = pysam.qualitystring_to_array(qual[::-1])
    else:
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array(qual)

    tags = [
        ("AS", int(map_data[off])),
        ("MD", str(map_data[off + 1])),
        ("ST", int(map_data[off + 2])),
    ]
    for i, tag in enumerate(["Yf", "Zf", "Yc", "Zc", "NS", "NC"]):
        tags.append((tag, int(map_data[off + 3 + i])))
    a.tags = tags
    return a


def create_unmapped_record(header, name, seq, qual, flag):
    a = pysam.AlignedSegment(header=header)
    a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = (
        name,
        flag,
        -1,
        -1,
        0,
    )
    a.query_sequence, a.query_qualities = seq, pysam.qualitystring_to_array(qual)
    return a


def _setup_output_bams(
    output_files, unmap_file, unified_header, ref_headers, stack, threads=1
):
    if len(output_files) > 1:
        bam_outs = {
            i: stack.enter_context(
                pysam.AlignmentFile(f, "wb", header=ref_headers[i], threads=threads)
            )
            for i, f in enumerate(output_files)
        }
        default_bam = bam_outs[len(output_files) - 1]
    else:
        default_bam = stack.enter_context(
            pysam.AlignmentFile(
                output_files[0], "wb", header=unified_header, threads=threads
            )
        )
        bam_outs = None
    bam_unmap = (
        stack.enter_context(
            pysam.AlignmentFile(
                unmap_file, "wb", header=unified_header, threads=threads
            )
        )
        if unmap_file
        else default_bam
    )
    return bam_outs, bam_unmap, default_bam


def map_file(
    r1_file,
    r2_file,
    ref_files,
    output_files,
    unmap_file=None,
    forward_library=True,
    max_mismatches=10,
    threads=8,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
    index_dir=None,
    index_only=False,
    batch_size=2000,
    orientation_filter=None,
):
    global _ORIENTATION_FILTER, _FORWARD_LIBRARY
    _ORIENTATION_FILTER, _FORWARD_LIBRARY = orientation_filter, forward_library
    ref_indices = [
        (None, None, str(i + 1) if len(ref_files) > 1 else "")
        for i in range(len(ref_files))
    ]
    i_dirs = (
        [index_dir]
        if index_dir and not isinstance(index_dir, (list, tuple))
        else (index_dir or [os.path.dirname(os.path.abspath(f)) for f in ref_files])
    )
    ref_indices = _build_and_check_indices(ref_files, ref_indices, i_dirs, index_only)
    if index_only:
        return
    mp_ctx = (
        mp.get_context("fork") if sys.platform != "win32" else mp.get_context("spawn")
    )
    stack = ExitStack()
    with stack:
        ref_headers = []
        all_sq = []
        seen_sq = set()
        for o_fa, _, _ in ref_indices:
            with pysam.FastaFile(o_fa) as fa:
                h = {
                    "HD": {"VN": "1.6", "SO": "unsorted"},
                    "SQ": [
                        {"LN": fa.get_reference_length(n), "SN": n}
                        for n in fa.references
                    ],
                }
                ref_headers.append(h)
                [
                    (all_sq.append(sq), seen_sq.add(sq["SN"]))
                    for sq in h["SQ"]
                    if sq["SN"] not in seen_sq
                ]
        unified_h = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": all_sq}
        global_rid_map = {sq["SN"]: i for i, sq in enumerate(all_sq)}
        global_rid_map["*"] = -1
        bam_outs, bam_unmap, def_bam = _setup_output_bams(
            output_files, unmap_file, unified_h, ref_headers, stack, threads=2
        )
        p_tot = p_map = 0
        s_time = time.time()
        with ProcessPoolExecutor(
            max_workers=max(1, threads),
            mp_context=mp_ctx,
            initializer=_init_worker,
            initargs=(ref_indices, orientation_filter, forward_library),
        ) as ex:
            futures = deque()
            batches = deque()
            it_reads = (
                read_paired_fastx(r1_file, r2_file) if r2_file else fastx_read(r1_file)
            )
            batch = []
            with Progress(
                SpinnerColumn(),
                TextColumn("[bold blue]Map[/bold blue]"),
                TextColumn("{task.fields[mapped]} / {task.fields[total_reads]} reads"),
                TextColumn("({task.fields[elapsed]})"),
                transient=False,
            ) as prg:
                task = prg.add_task(
                    "Mapping", total=None, mapped="0", total_reads="0", elapsed="0.00s"
                )

                def flush_done():
                    nonlocal p_tot, p_map
                    while futures and futures[0].done():
                        res = futures.popleft().result()
                        orig_batch = batches.popleft()
                        for i, (ridx, m_res) in enumerate(res):
                            p_tot += 1
                            info = orig_batch[i]
                            if ridx is not None and m_res:
                                p_map += 1
                                target = bam_outs[ridx] if bam_outs else def_bam
                                for j, item in enumerate(m_res):
                                    if r2_file:
                                        if item[1]:
                                            target.write(
                                                create_bam_record(
                                                    target.header,
                                                    item[1],
                                                    info[0][0],
                                                    info[0][1],
                                                    info[0][2],
                                                    j > 0,
                                                    global_rid_map,
                                                    1,
                                                )
                                            )
                                        if item[2]:
                                            target.write(
                                                create_bam_record(
                                                    target.header,
                                                    item[2],
                                                    info[1][0],
                                                    info[1][1],
                                                    info[1][2],
                                                    j > 0,
                                                    global_rid_map,
                                                    2,
                                                )
                                            )
                                    else:
                                        target.write(
                                            create_bam_record(
                                                target.header,
                                                item,
                                                info[0],
                                                info[1],
                                                info[2],
                                                j > 0,
                                                global_rid_map,
                                                0,
                                            )
                                        )
                            else:
                                if r2_file:
                                    bam_unmap.write(
                                        create_unmapped_record(
                                            bam_unmap.header,
                                            info[0][0],
                                            info[0][1],
                                            info[0][2],
                                            77,
                                        )
                                    )
                                    bam_unmap.write(
                                        create_unmapped_record(
                                            bam_unmap.header,
                                            info[1][0],
                                            info[1][1],
                                            info[1][2],
                                            141,
                                        )
                                    )
                                else:
                                    bam_unmap.write(
                                        create_unmapped_record(
                                            bam_unmap.header,
                                            info[0],
                                            info[1],
                                            info[2],
                                            4,
                                        )
                                    )
                        prg.update(
                            task,
                            mapped=f"{p_map:,}",
                            total_reads=f"{p_tot:,}",
                            elapsed=format_duration(time.time() - s_time),
                        )

                for reads in it_reads:
                    if r2_file:
                        if reads[0].name.split()[0].rstrip("/1").rstrip("/2") != reads[
                            1
                        ].name.split()[0].rstrip("/1").rstrip("/2"):
                            raise ValueError(
                                f"Order mismatch: {reads[0].name} vs {reads[1].name}"
                            )
                        batch.append(
                            (
                                (reads[0].name, reads[0].sequence, reads[0].quality),
                                (reads[1].name, reads[1].sequence, reads[1].quality),
                            )
                        )
                    else:
                        batch.append((reads.name, reads.sequence, reads.quality))
                    if len(batch) >= batch_size:
                        batches.append(list(batch))
                        futures.append(
                            ex.submit(
                                _map_batch_worker,
                                batches[-1],
                                r2_file is not None,
                                max_mismatches,
                                min_alignment_length,
                                min_mapping_ratio,
                            )
                        )
                        batch.clear()
                        while len(futures) >= threads * 2:
                            wait(futures, return_when=FIRST_COMPLETED)
                            flush_done()
                if batch:
                    batches.append(list(batch))
                    futures.append(
                        ex.submit(
                            _map_batch_worker,
                            batches[-1],
                            r2_file is not None,
                            max_mismatches,
                            min_alignment_length,
                            min_mapping_ratio,
                        )
                    )
                while futures:
                    wait(futures, return_when=FIRST_COMPLETED)
                    flush_done()
    print(f"\n✅ Mapping completed! Output saved to: {output_files}")
