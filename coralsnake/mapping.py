#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-08 20:32


import atexit
import multiprocessing as mp
import os
import random
import shutil
import tempfile
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack

import pysam
from bwamem import BwaAligner, BwaIndexer, fastx_read, read_paired_fastx
from rich.progress import Progress, SpinnerColumn, TextColumn

from . import seqops
from .utils import (
    convert_file_realtime,
    format_duration,
    km_conversion,
    mk_conversion,
    reverse_complement,
)


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

    def build_orig(): indexer1.build_index(orig_fa, prefix=orig_prefix)
    def build_mk(): indexer2.build_index(mk_fa, prefix=mk_prefix)
    def do_conversion(): convert_file_realtime(ref_file, mk_fa, "AC", "GT")

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
            if fut_orig.done() and fut_mk.done(): break
            time.sleep(poll_interval)
        fut_orig.result(); fut_mk.result()

    on_update("Done", 100.0, 100.0)
    return orig_fa, mk_prefix


def _map_batch_worker(batch, paired, max_mismatches, min_alignment_length, min_mapping_ratio):
    results = []
    if not paired:
        seqs = [item[1] for item in batch]
        conv1 = seqops.batch_base_conversion(seqs, "AC", "GT") if _FORWARD_LIBRARY else seqops.batch_base_conversion(seqs, "GT", "AC")
        conv2 = seqops.batch_base_conversion(seqs, "GT", "AC") if _FORWARD_LIBRARY else seqops.batch_base_conversion(seqs, "AC", "GT")
        for i, (name1, seq1, qua1) in enumerate(batch):
            ref_idx, mapping_result = run_mapping_se(name1, seq1, qua1, conv1[i], conv2[i], _FORWARD_LIBRARY, max_mismatches, min_alignment_length, min_mapping_ratio)
            results.append(((name1, seq1, qua1), (ref_idx, mapping_result)))
    else:
        seqs1 = [item[0][1] for item in batch]; seqs2 = [item[1][1] for item in batch]
        rc_s1 = [seqops.reverse_complement(s) for s in seqs1]; rc_s2 = [seqops.reverse_complement(s) for s in seqs2]
        # Pre-convert for native orientation only (most common case)
        # Orientation 1: Read1 converted, Read2-RC converted
        # Orientation 2: Read1-RC converted, Read2 converted
        if _FORWARD_LIBRARY:
            c1_r1 = seqops.batch_base_conversion(seqs1, "AC", "GT")
            c1_r2 = seqops.batch_base_conversion(rc_s2, "AC", "GT")
            c2_r1 = seqops.batch_base_conversion(rc_s1, "AC", "GT")
            c2_r2 = seqops.batch_base_conversion(seqs2, "AC", "GT")
        else:
            c1_r1 = seqops.batch_base_conversion(rc_s1, "AC", "GT")
            c1_r2 = seqops.batch_base_conversion(seqs2, "AC", "GT")
            c2_r1 = seqops.batch_base_conversion(seqs1, "AC", "GT")
            c2_r2 = seqops.batch_base_conversion(rc_s2, "AC", "GT")
        conv_data = (c1_r1, c1_r2, c2_r1, c2_r2)
        for i, ((name1, seq1, qua1), (name2, seq2, qua2)) in enumerate(batch):
            ref_idx, mapping_result = run_mapping_pe(name1, seq1, seq2, qua1, qua2, i, conv_data, rc_s1[i], rc_s2[i], _FORWARD_LIBRARY, max_mismatches, min_alignment_length, min_mapping_ratio)
            results.append(((name1, seq1, qua1, name2, seq2, qua2), (ref_idx, mapping_result)))
    return results


def _check_reference_length(aligner, min_length=100):
    for i in range(aligner.index.bns.n_seqs):
        if aligner.index.bns.anns[i].len >= min_length: return True
    return False


def run_mapping_se(name, seq1, qua1, s1_c1, s1_c2, forward_library, max_mismatches, min_alignment_length, min_mapping_ratio):
    rc_seq1 = rc_qua1 = None
    for ref_idx in range(len(_ALIGNERS_MK)):
        if not _check_reference_length(_ALIGNERS_MK[ref_idx]): continue
        mapped = []
        oris = [1, 2] if _ORIENTATION_FILTER is None else [_ORIENTATION_FILTER]
        for orientation in oris:
            s1_c = s1_c1 if orientation == 1 else s1_c2
            hits = _ALIGNERS_MK[ref_idx].align(s1_c, min_mapq=1)
            for h in hits:
                if (h[2]-h[1]) < min_alignment_length: continue
                ref = _ALIGNERS_ORIG[ref_idx].seq(h[0], h[1], h[2]+100)
                r1_rev = h[3] == -1
                if r1_rev:
                    if rc_seq1 is None: rc_seq1, rc_qua1 = seqops.reverse_complement(seq1), qua1[::-1]
                    s, q, flag = rc_seq1, rc_qua1, 16
                else: s, q, flag = seq1, qua1, 0
                res = seqops.score_and_tag(h[7], s, ref, (orientation == 1) ^ r1_rev)
                if res[1] > max_mismatches // 2 or (len(seq1) * min_mapping_ratio > (len(seq1) - res[1])): continue
                m1 = [name, flag, h[0], h[1]+1, score_to_mapq(res[0]), h[7], "*", 0, 0, s, q, ("MD", res[2]), ("ST", orientation), ("AS", res[0]), ("Yf", res[3]), ("Zf", res[4]), ("Yc", res[5]), ("Zc", res[6]), ("NS", res[7]), ("NC", res[8])]
                mapped.append([res[0], m1])
        if mapped:
            random.shuffle(mapped); mapped = sorted(mapped, key=lambda x: x[0], reverse=True)
            return (ref_idx, mapped)
    return (None, [])


def run_mapping_pe(name, seq1, seq2, qua1, qua2, idx, conv_data, rc_s1, rc_s2, forward_library, max_mismatches, min_alignment_length, min_mapping_ratio):
    q1_l, q2_l = len(seq1), len(seq2)
    rc_q1 = rc_q2 = None
    c1_r1, c1_r2, c2_r1, c2_r2 = conv_data
    all_res = []
    for ref_idx in range(len(_ALIGNERS_MK)):
        if not _check_reference_length(_ALIGNERS_MK[ref_idx]): continue
        mapped = []
        oris = [1, 2] if _ORIENTATION_FILTER is None else [_ORIENTATION_FILTER]
        for orientation in oris:
            is_o1 = orientation == 1
            s1_c, s2_c = (c1_r1[idx], c1_r2[idx]) if is_o1 else (c2_r1[idx], c2_r2[idx])
            hits = _ALIGNERS_MK[ref_idx].align(s1_c, s2_c, min_mapq=1)
            for h1, h2, is_p, isize in hits:
                res1 = res2 = None
                if h1 and (h1[2]-h1[1]) >= min_alignment_length:
                    ref1 = _ALIGNERS_ORIG[ref_idx].seq(h1[0], h1[1], h1[2]+100)
                    r1_rev = (h1[3]==-1 if forward_library else h1[3]==1) if is_o1 else (h1[3]==1 if forward_library else h1[3]==-1)
                    if r1_rev:
                        if rc_q1 is None: rc_q1 = qua1[::-1]
                        s1, q1, f1 = rc_s1, rc_q1, 16
                    else: s1, q1, f1 = seq1, qua1, 0
                    res1 = seqops.score_and_tag(h1[7], s1, ref1, (is_o1 if forward_library else (not is_o1)) ^ r1_rev)
                if h2 and (h2[2]-h2[1]) >= min_alignment_length:
                    ref2 = _ALIGNERS_ORIG[ref_idx].seq(h2[0], h2[1], h2[2]+100)
                    r2_rev = (h2[3]==1 if forward_library else h2[3]==-1) if is_o1 else (h2[3]==-1 if forward_library else h2[3]==1)
                    if r2_rev:
                        if rc_q2 is None: rc_q2 = qua2[::-1]
                        s2, q2, f2 = rc_s2, rc_q2, 16
                    else: s2, q2, f2 = seq2, qua2, 0
                    res2 = seqops.score_and_tag(h2[7], s2, ref2, (not (is_o1 if forward_library else (not is_o1))) ^ r2_rev)
                if res1 or res2:
                    if ( (res1[1] if res1 else 0) + (res2[1] if res2 else 0) ) > max_mismatches: continue
                    if res1 and (q1_l * min_mapping_ratio > (q1_l - res1[1])): continue
                    if res2 and (q2_l * min_mapping_ratio > (q2_l - res2[1])): continue
                    mq = score_to_mapq(min(res1[0] if res1 else 99, res2[0] if res2 else 99))
                    m1 = m2 = None
                    if res1:
                        t1 = [("ST", orientation), ("MD", res1[2]), ("AS", res1[0]), ("Yf", res1[3]), ("Zf", res1[4]), ("Yc", res1[5]), ("Zc", res1[6]), ("NS", res1[7]), ("NC", res1[8])]
                        if h2: m1 = [name, 67+f1+(32 if h2[3]==-1 else 0)+(2 if is_p else 0), h1[0], h1[1]+1, mq, h1[7], h2[0], h2[1]+1, (isize if h1[1]<=h2[1] else -isize), s1, q1] + t1
                        else: m1 = [name, 73+f1, h1[0], h1[1]+1, mq, h1[7], "*", 0, 0, s1, q1] + t1
                    if res2:
                        t2 = [("ST", orientation), ("MD", res2[2]), ("AS", res2[0]), ("Yf", res2[3]), ("Zf", res2[4]), ("Yc", res2[5]), ("Zc", res2[6]), ("NS", res2[7]), ("NC", res2[8])]
                        if h1: m2 = [name, 131+f2+(32 if h1[3]==-1 else 0)+(2 if is_p else 0), h2[0], h2[1]+1, mq, h2[7], h1[0], h1[1]+1, (isize if h2[1]<=h1[1] else -isize), s2, q2] + t2
                        else: m2 = [name, 137+f2, h2[0], h2[1]+1, mq, h2[7], "*", 0, 0, s2, q2] + t2
                    mapped.append([(res1[0] if res1 else 0)+(res2[0] if res2 else 0), m1, m2])
        if mapped:
            random.shuffle(mapped); mapped = sorted(mapped, key=lambda x: x[0], reverse=True)
            if any(r[1] and r[2] for r in mapped): return (ref_idx, mapped)
            all_res.append((ref_idx, mapped))
        else: all_res.append((ref_idx, []))
    for r_idx, rs in all_res:
        if rs: return (r_idx, rs)
    return (None, [])


_ALIGNERS_ORIG = []; _ALIGNERS_MK = []; _ORIENTATION_FILTER = None; _FORWARD_LIBRARY = None

def _init_worker(ref_indices, orientation_filter, forward_library):
    global _ALIGNERS_ORIG, _ALIGNERS_MK, _ORIENTATION_FILTER, _FORWARD_LIBRARY
    _ORIENTATION_FILTER, _FORWARD_LIBRARY = orientation_filter, forward_library
    _ALIGNERS_ORIG = []; _ALIGNERS_MK = []
    for o_fa, mk_pre, _ in ref_indices:
        opts = {"min_seed_len": 14, "max_occ": 1000, "softclip_supplementary": True, "mark_secondary": True, "clip_penalties": (6,6), "unpaired_penalty": 24, "min_score": 20, "insert_model": (80,60,450)}
        _ALIGNERS_ORIG.append(BwaAligner(os.path.splitext(o_fa)[0], **opts))
        _ALIGNERS_MK.append(BwaAligner(mk_pre, **opts))


def _build_and_check_indices(ref_files, ref_indices, index_dirs, index_only):
    with Progress(SpinnerColumn(), TextColumn("[bold green]Index[/bold green]"), TextColumn("{task.description}"), TextColumn("| [yellow]{task.fields[status]}[/yellow] | [cyan]{task.fields[progress]}[/cyan]"), transient=False) as progress:
        task = progress.add_task(f"Indices...", total=None, status="Checking", progress="0%")
        for i in range(len(ref_indices)):
            ref_suffix = ref_indices[i][2]
            d = index_dirs[i] if len(index_dirs) == len(ref_indices) else index_dirs[0]
            o_fa, mk_pre = os.path.join(d, f"ref{ref_suffix}.orig.fa"), os.path.join(d, f"ref{ref_suffix}.mk")
            if index_only or not (os.path.exists(o_fa) and all(os.path.exists(mk_pre+e) for e in [".amb",".ann",".bwt",".pac",".sa"])):
                def on_up(s,p1,p2): progress.update(task, status=f"Ref {i+1}: {s}", progress=f"{(p1+p2)/2:.1f}%")
                idx0, idx_mk = _build_indices_with_progress(ref_files[i] if i < len(ref_files) else None, d, on_up, ref_suffix=ref_suffix)
                ref_indices[i] = (idx0, idx_mk, ref_suffix)
            else: ref_indices[i] = (o_fa, mk_pre, ref_suffix)
        progress.update(task, description="✓ Indices ready", status="Done", progress="100%")
        return ref_indices


def _process_and_map_reads(r1_file, r2_file, ref_indices, write_func, batch_size, threads, max_mismatches, min_alignment_length, min_mapping_ratio, orientation_filter, forward_library):
    paired = r2_file is not None
    it_reads = read_paired_fastx(r1_file, r2_file) if paired else fastx_read(r1_file)
    batch = []
    with ProcessPoolExecutor(max_workers=max(1, threads), mp_context=mp.get_context("spawn"), initializer=_init_worker, initargs=(ref_indices, orientation_filter, forward_library)) as ex:
        futures = deque()
        for reads in it_reads:
            if paired:
                if reads[0].name.split()[0].rstrip("/1").rstrip("/2") != reads[1].name.split()[0].rstrip("/1").rstrip("/2"):
                    raise ValueError(f"Order mismatch: {reads[0].name} vs {reads[1].name}")
                batch.append(((reads[0].name, reads[0].sequence, reads[0].quality), (reads[1].name, reads[1].sequence, reads[1].quality)))
            else: batch.append((reads.name, reads.sequence, reads.quality))
            if len(batch) >= batch_size:
                futures.append(ex.submit(_map_batch_worker, list(batch), paired, max_mismatches, min_alignment_length, min_mapping_ratio))
                batch.clear()
            while futures and futures[0].done(): write_func(futures.popleft().result())
        if batch: futures.append(ex.submit(_map_batch_worker, list(batch), paired, max_mismatches, min_alignment_length, min_mapping_ratio))
        while futures: write_func(futures.popleft().result())


def _write_mapped_record(paired, item, read_info, target_bam, i):
    if paired and len(item) == 3:
        if item[1]: target_bam.write(create_bam_record(target_bam.header, item[1], i > 0))
        if item[2]: target_bam.write(create_bam_record(target_bam.header, item[2], i > 0))
    else:
        target_bam.write(create_bam_record(target_bam.header, item[1], i > 0))


def _write_unmapped_record(paired, read_info, bam_unmap):
    if paired:
        bam_unmap.write(create_unmapped_record(bam_unmap.header, read_info[0], read_info[1], read_info[2], 77))
        bam_unmap.write(create_unmapped_record(bam_unmap.header, read_info[3], read_info[4], read_info[5], 141))
    else: bam_unmap.write(create_unmapped_record(bam_unmap.header, read_info[0], read_info[1], read_info[2], 4))


def _setup_output_bams(output_files, unmap_file, unified_header, ref_headers, stack):
    if len(output_files) > 1:
        bam_outs = {i: stack.enter_context(pysam.AlignmentFile(f, "wb", header=ref_headers[i])) for i, f in enumerate(output_files)}
        default_bam = bam_outs[len(output_files)-1]
    else:
        default_bam = stack.enter_context(pysam.AlignmentFile(output_files[0], "wb", header=unified_h)) if 'unified_h' in locals() else stack.enter_context(pysam.AlignmentFile(output_files[0], "wb", header=unified_header))
        bam_outs = None
    bam_unmap = stack.enter_context(pysam.AlignmentFile(unmap_file, "wb", header=locals().get('unified_h', locals().get('unified_header')))) if unmap_file else default_bam
    return bam_outs, bam_unmap, default_bam


def score_to_mapq(score): return max(0, min(60, score))

def create_bam_record(header, map_data, is_secondary):
    a = pysam.AlignedSegment(header=header)
    a.query_name, a.flag = map_data[0], map_data[1] + (256 if is_secondary else 0)
    a.reference_name, a.reference_start, a.mapping_quality, a.cigarstring = map_data[2], map_data[3]-1, map_data[4], map_data[5]
    a.next_reference_name = map_data[6]
    a.next_reference_start = map_data[7]-1 if map_data[7]>0 else 0
    a.template_length, a.query_sequence, a.query_qualities = map_data[8], map_data[9], pysam.qualitystring_to_array(map_data[10])
    for t, v in map_data[11:]: a.set_tag(t, v)
    return a

def create_unmapped_record(header, name, seq, qual, flag):
    a = pysam.AlignedSegment(header=header)
    a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = name, flag, -1, -1, 0
    a.cigarstring, a.next_reference_id, a.next_reference_start, a.template_length = None, -1, -1, 0
    a.query_sequence, a.query_qualities = seq, pysam.qualitystring_to_array(qual)
    return a


def map_file(r1_file, r2_file, ref_files, output_files, unmap_file=None, forward_library=True, max_mismatches=10, threads=8, min_alignment_length=20, min_mapping_ratio=0.5, index_dir=None, index_only=False, batch_size=5000, orientation_filter=None):
    global _ORIENTATION_FILTER, _FORWARD_LIBRARY
    _ORIENTATION_FILTER, _FORWARD_LIBRARY = orientation_filter, forward_library
    ref_indices = [(None, None, str(i+1) if len(ref_files)>1 else "") for i in range(len(ref_files))]
    i_dirs = [index_dir] if index_dir and not isinstance(index_dir, (list,tuple)) else (index_dir or [os.path.dirname(os.path.abspath(f)) for f in ref_files])
    ref_indices = _build_and_check_indices(ref_files, ref_indices, i_dirs, index_only)
    if index_only: return
    stack = ExitStack()
    with stack:
        ref_headers = []; all_sq = []; seen_sq = set()
        for o_fa, _, _ in ref_indices:
            with pysam.FastaFile(o_fa) as fa:
                h = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"LN": fa.get_reference_length(n), "SN": n} for n in fa.references]}
                ref_headers.append(h)
                for sq in h["SQ"]:
                    if sq["SN"] not in seen_sq: all_sq.append(sq); seen_sq.add(sq["SN"])
        unified_h = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": all_sq}
        bam_outs, bam_unmap, def_bam = _setup_output_bams(output_files, unmap_file, unified_h, ref_headers, stack)
        p_tot = p_map = 0; s_time = time.time()
        def wr(res):
            nonlocal p_tot, p_map
            for info, (ridx, m_res) in res:
                p_tot += 1
                if ridx is not None and m_res:
                    p_map += 1; target = bam_outs[ridx] if bam_outs else def_bam
                    for i, item in enumerate(m_res): _write_mapped_record(r2_file is not None, item, info, target, i)
                else: _write_unmapped_record(r2_file is not None, info, bam_unmap)
        with Progress(SpinnerColumn(), TextColumn("[bold blue]Map[/bold blue]"), TextColumn("{task.fields[mapped]} / {task.fields[total_reads]} reads"), TextColumn("({task.fields[elapsed]})"), transient=False) as prg:
            task = prg.add_task("Mapping", total=None, mapped="0", total_reads="0", elapsed="0.00s")
            def pw(res):
                wr(res); prg.update(task, mapped=f"{p_map:,}", total_reads=f"{p_tot:,}", elapsed=format_duration(time.time()-s_time))
            _process_and_map_reads(r1_file, r2_file, ref_indices, pw, batch_size, threads, max_mismatches, min_alignment_length, min_mapping_ratio, orientation_filter, forward_library)
    print(f"\n✅ Mapping completed! Output saved to: {output_files}")
