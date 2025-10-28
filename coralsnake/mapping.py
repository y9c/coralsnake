#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack

import pysam
from bwamem import BwaAligner, BwaIndexer, fastx_read, read_paired_fastx
from rich.progress import Progress, SpinnerColumn, TextColumn

from . import seqops
from .utils import convert_file_realtime, format_duration, km_conversion, mk_conversion


## (removed) async version: asyncio.run() blocks Rich progress updates


def _build_indices_with_progress(
    ref_file: str,
    index_base_dir: str,
    on_update,
    poll_interval: float = 0.1,
):
    """Build indices with progress using threads (not async - Rich doesn't play well with asyncio.run)."""
    from concurrent.futures import ThreadPoolExecutor
    
    os.makedirs(index_base_dir, exist_ok=True)
    # ORIG
    orig_fa = os.path.join(index_base_dir, "ref.orig.fa")
    if os.path.abspath(ref_file) != os.path.abspath(orig_fa):
        shutil.copyfile(ref_file, orig_fa)
    orig_prefix = os.path.splitext(orig_fa)[0]
    # MK
    mk_fa = os.path.join(index_base_dir, "ref.mk.fa")
    mk_prefix = os.path.splitext(mk_fa)[0]

    # Start ORIG index in parallel with conversion (they're independent)
    import time
    
    indexer1 = BwaIndexer()
    indexer2 = BwaIndexer()

    def build_orig():
        indexer1.build_index(orig_fa, prefix=orig_prefix, capture_progress=True)

    def build_mk():
        indexer2.build_index(mk_fa, prefix=mk_prefix, capture_progress=True)
    
    def do_conversion():
        convert_file_realtime(ref_file, mk_fa, "AC", "GT")

    with ThreadPoolExecutor(max_workers=3) as ex:
        # Start ORIG index and conversion in parallel
        fut_orig = ex.submit(build_orig)
        fut_conv = ex.submit(do_conversion)
        
        # Poll while conversion runs
        conv_status = "Converting..."
        while not fut_conv.done():
            # If ORIG is done during conversion, mark it complete
            p1 = 100.0 if fut_orig.done() else (indexer1.progress_percent or 0.0)
            on_update(conv_status, p1, 0.0)
            time.sleep(poll_interval)
        
        fut_conv.result()
        conv_status = "Converted"
        
        # Start MK index after conversion completes
        fut_mk = ex.submit(build_mk)
        
        # Poll both indices
        while not (fut_orig.done() and fut_mk.done()):
            # If future is done but progress still 0, it completed too fast - set to 100%
            p1 = 100.0 if fut_orig.done() else (indexer1.progress_percent or 0.0)
            p2 = 100.0 if fut_mk.done() else (indexer2.progress_percent or 0.0)
            on_update(conv_status, p1, p2)
            time.sleep(poll_interval)
        
        fut_orig.result()
        fut_mk.result()

    on_update("Done", 100.0, 100.0)
    return orig_fa, mk_prefix


## (removed) _ensure_indices: inlined async index building with progress in map_file


def _map_batch_worker(
    batch,
    paired,
    fwd_lib,
    max_mismatches,
    min_alignment_length,
    min_mapping_ratio,
):
    """Map a batch of reads; return one run_mapping result per input read."""
    # Use per-process cached resources initialized by _init_worker

    # timing removed
    results = []
    if not paired:
        for name1, seq1, qua1 in batch:
            results.append(
                run_mapping_se(
                    name1,
                    seq1,
                    qua1,
                    fwd_lib,
                    max_mismatches,
                    min_alignment_length,
                    min_mapping_ratio,
                )
            )
    else:
        for (name1, seq1, qua1), (name2, seq2, qua2) in batch:
            base1 = name1.split()[0].rstrip("/1").rstrip("/2")
            base2 = name2.split()[0].rstrip("/1").rstrip("/2")
            if base1 != base2:
                raise ValueError(f"r1 and r2 not in the same order: {name1} vs {name2}")
            results.append(
                run_mapping_pe(
                    name1,
                    seq1,
                    seq2,
                    qua1,
                    qua2,
                    fwd_lib,
                    max_mismatches,
                    min_alignment_length,
                    min_mapping_ratio,
                )
            )

    return results


# Per-process cached aligners (initialized once per worker)
_ALIGNER_ORIG = None
_ALIGNER_MK = None


def _init_worker(orig_fa, mk_index_prefix, threads):
    global _ALIGNER_ORIG, _ALIGNER_MK
    _ALIGNER_ORIG = BwaAligner(
        os.path.splitext(orig_fa)[0],
        softclip_supplementary=True,
        mark_secondary=True,
        clip_penalties=(6, 6),
        unpaired_penalty=24,
        min_score=20,
        insert_model=(80, 60, 450),
    )
    _ALIGNER_MK = BwaAligner(
        mk_index_prefix,
        softclip_supplementary=True,
        mark_secondary=True,
        clip_penalties=(6, 6),
        unpaired_penalty=24,
        min_score=20,
        insert_model=(80, 60, 450),
    )


def _revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def find_properly_paired_hits(hits, fwd=True):
    """Find read1/read2 hit pairs on the same contig, opposite strands, within 1 kb."""
    parsed_hits = []
    # group by ref_name and separate read 1 and read 2
    ref_name_hits = {}
    for hit in hits:
        if hit.ctg not in ref_name_hits:
            ref_name_hits[hit.ctg] = [[], []]
        ref_name_hits[hit.ctg][hit.read_num - 1].append(hit)
    for hits in ref_name_hits.values():
        if len(hits[0]) > 0 and len(hits[1]) > 0:
            for hit1 in hits[0]:
                for hit2 in hits[1]:
                    if hit1.strand + hit2.strand == 0:
                        if fwd:
                            if hit1.r_st < hit2.r_en and hit2.r_en - hit1.r_st < 1000:
                                parsed_hits.append((hit1, hit2))
                        else:
                            if hit1.r_en > hit2.r_st and hit1.r_en - hit2.r_st < 1000:
                                parsed_hits.append((hit1, hit2))

    return parsed_hits


def cal_md_and_tag(cigar, seq, ref, fwd):
    """Compute MD tag and conversion stats: return (md, yf, zf, yc, zc, ns, nc)."""
    # Use optimized C implementation
    return seqops.cal_md_and_tag(cigar, seq, ref, fwd)


def calculate_directional_score(cigar, seq, ref, is_orientation1):
    """Score alignment with conversion awareness: return (score, wrong_conversions, bad_mismatches)."""
    # Use optimized C implementation
    return seqops.calculate_directional_score(cigar, seq, ref, is_orientation1)


def filter_hits(hits, seq1, seq2, min_alignment_length=20, min_mapping_ratio=0.5):
    """Keep hits with mapq>0, blen>=min_alignment_length, mlen>=min_mapping_ratio*qlen."""
    filtered_hits = []
    for hit in hits:
        q_len = len(seq1) if hit.read_num == 1 else len(seq2)
        if (
            hit.mapq > 0
            and hit.blen > min_alignment_length
            and hit.mlen > min_mapping_ratio * q_len
        ):
            filtered_hits.append(hit)
    return filtered_hits


def run_mapping_se(
    name,
    seq1,
    qua1,
    fwd_lib=True,
    max_mismatches=10,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
):
    """Map one single-end read and return scored alignments."""
    mapped = []
    for orientation in [1, 2]:
        # Build converted read
        if orientation == 1:
            seq1_conv = mk_conversion(seq1) if fwd_lib else km_conversion(seq1)
        else:
            seq1_conv = km_conversion(seq1) if fwd_lib else mk_conversion(seq1)

        # Iterate hits
        # Align converted read to MK reference using BWA
        hits = tuple(_ALIGNER_MK.align(seq1_conv))
        for h in hits:
            try:
                h.read_num = 1
            except Exception:
                pass
        for hit in filter_hits(
            hits,
            seq1,
            None,
            min_alignment_length,
            min_mapping_ratio,
        ):
            ref = _ALIGNER_ORIG.seq(hit.ctg, hit.r_st, hit.r_en)
            read_reverse = hit.strand == -1
            if read_reverse:
                flag = 16
                s = _revcomp(seq1)
                q = qua1[::-1]
            else:
                flag = 0
                s = seq1
                q = qua1

            cigar_str = hit.cigar_str
            cigar = hit.cigar
            if hit.q_st > 0:
                cigar_str = f"{hit.q_st}S" + cigar_str
                cigar = [[hit.q_st, 4]] + cigar
            if hit.q_en < len(s):
                cigar_str = cigar_str + f"{len(s) - hit.q_en}S"
                cigar = cigar + [[len(s) - hit.q_en, 4]]

            is_orientation1 = orientation == 1
            score, _wrong_conv, bad_mm = calculate_directional_score(
                cigar, s, ref, is_orientation1
            )
            if bad_mm > max_mismatches // 2:
                continue

            md, yf, zf, yc, zc, ns, nc = cal_md_and_tag(cigar, s, ref, is_orientation1)
            mapq = min(60, score)
            tags = [
                ("MD", md),
                ("ST", orientation),
                ("AS", score),
                ("Yf", yf),
                ("Zf", zf),
                ("Yc", yc),
                ("Zc", zc),
                ("NS", ns),
                ("NC", nc),
            ]
            map1 = [
                name,
                flag,
                hit.ctg,
                hit.r_st + 1,
                mapq,
                cigar_str,
                "*",
                0,
                0,
                s,
                q,
            ] + tags
            mapped.append([score, map1])

    random.shuffle(mapped)
    mapped = sorted(mapped, key=lambda x: x[0], reverse=True)
    return mapped


def run_mapping_pe(
    name,
    seq1,
    seq2,
    qua1,
    qua2,
    fwd_lib=True,
    max_mismatches=10,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
):
    """Map one paired-end read and return scored alignment pairs."""
    mapped = []
    for orientation in [1, 2]:
        # Build converted reads
        if orientation == 1:
            if fwd_lib:
                seq1_conv = mk_conversion(seq1)
                seq2_conv = km_conversion(seq2)
            else:
                seq1_conv = km_conversion(seq1)
                seq2_conv = mk_conversion(seq2)
        else:
            if fwd_lib:
                seq1_conv = km_conversion(seq1)
                seq2_conv = mk_conversion(seq2)
            else:
                seq1_conv = mk_conversion(seq1)
                seq2_conv = km_conversion(seq2)

        # Align both converted reads independently to MK reference
        hits1 = tuple(_ALIGNER_MK.align(seq1_conv))
        for h in hits1:
            try:
                h.read_num = 1
            except Exception:
                pass
        hits2 = tuple(_ALIGNER_MK.align(seq2_conv))
        for h in hits2:
            try:
                h.read_num = 2
            except Exception:
                pass
        combined_hits = hits1 + hits2
        for hit1, hit2 in find_properly_paired_hits(
            filter_hits(
                combined_hits,
                seq1,
                seq2,
                min_alignment_length,
                min_mapping_ratio,
            ),
            fwd=True,
        ):
            tlen = max(hit1.r_en, hit2.r_en) - min(hit1.r_st, hit2.r_st)
            ref1 = _ALIGNER_ORIG.seq(hit1.ctg, hit1.r_st, hit1.r_en)
            ref2 = _ALIGNER_ORIG.seq(hit2.ctg, hit2.r_st, hit2.r_en)
            read1_reverse = hit1.strand == -1
            read2_reverse = hit2.strand == -1
            if read1_reverse:
                s1 = _revcomp(seq1)
                q1 = qua1[::-1]
            else:
                s1 = seq1
                q1 = qua1
            if read2_reverse:
                s2 = _revcomp(seq2)
                q2 = qua2[::-1]
            else:
                s2 = seq2
                q2 = qua2

            if read1_reverse and not read2_reverse:
                flag1, flag2 = 83, 163
            elif not read1_reverse and read2_reverse:
                flag1, flag2 = 99, 147
            else:
                flag1, flag2 = 67, 131

            cigar_str1 = hit1.cigar_str
            cigar1 = hit1.cigar
            if hit1.q_st > 0:
                cigar_str1 = f"{hit1.q_st}S" + cigar_str1
                cigar1 = [[hit1.q_st, 4]] + cigar1
            if hit1.q_en < len(s1):
                cigar_str1 = cigar_str1 + f"{len(s1) - hit1.q_en}S"
                cigar1 = cigar1 + [[len(s1) - hit1.q_en, 4]]
            cigar_str2 = hit2.cigar_str
            cigar2 = hit2.cigar
            if hit2.q_st > 0:
                cigar_str2 = f"{hit2.q_st}S" + cigar_str2
                cigar2 = [[hit2.q_st, 4]] + cigar2
            if hit2.q_en < len(s2):
                cigar_str2 = cigar_str2 + f"{len(s2) - hit2.q_en}S"
                cigar2 = cigar2 + [[len(s2) - hit2.q_en, 4]]

            is_orientation1 = orientation == 1
            score1, _w1, bad_mm1 = calculate_directional_score(
                cigar1, s1, ref1, is_orientation1
            )
            score2, _w2, bad_mm2 = calculate_directional_score(
                cigar2, s2, ref2, is_orientation1
            )
            if (bad_mm1 + bad_mm2) > max_mismatches:
                continue

            md1, yf1, zf1, yc1, zc1, ns1, nc1 = cal_md_and_tag(
                cigar1, s1, ref1, is_orientation1
            )
            md2, yf2, zf2, yc2, zc2, ns2, nc2 = cal_md_and_tag(
                cigar2, s2, ref2, is_orientation1
            )
            combined_score = score1 + score2
            mapq = min(60, min(score1, score2))
            common_tags = [("ST", orientation)]
            tags1 = common_tags + [
                ("MD", md1),
                ("AS", score1),
                ("Yf", yf1),
                ("Zf", zf1),
                ("Yc", yc1),
                ("Zc", zc1),
                ("NS", ns1),
                ("NC", nc1),
            ]
            map1 = [
                name,
                flag1,
                hit1.ctg,
                hit1.r_st + 1,
                mapq,
                cigar_str1,
                hit2.ctg,
                hit2.r_st + 1,
                tlen,
                s1,
                q1,
            ] + tags1
            tags2 = common_tags + [
                ("MD", md2),
                ("AS", score2),
                ("Yf", yf2),
                ("Zf", zf2),
                ("Yc", yc2),
                ("Zc", zc2),
                ("NS", ns2),
                ("NC", nc2),
            ]
            map2 = [
                name,
                flag2,
                hit2.ctg,
                hit2.r_st + 1,
                mapq,
                cigar_str2,
                hit1.ctg,
                hit1.r_st + 1,
                -tlen,
                s2,
                q2,
            ] + tags2
            mapped.append([combined_score, map1, map2])

    random.shuffle(mapped)
    mapped = sorted(mapped, key=lambda x: x[0], reverse=True)
    return mapped


## (removed) run_mapping: worker dispatches directly to SE/PE


def create_bam_record(header, map_data, is_secondary):
    """Create a pysam.AlignedSegment from mapping data (SAM fields + tags)."""
    a = pysam.AlignedSegment(header=header)
    a.query_name = map_data[0]
    a.flag = map_data[1] + (256 if is_secondary else 0)
    a.reference_name = map_data[2]
    a.reference_start = map_data[3] - 1
    a.mapping_quality = map_data[4]
    a.cigarstring = map_data[5]
    a.next_reference_name = map_data[6]
    a.next_reference_start = map_data[7] - 1 if map_data[7] > 0 else 0
    a.template_length = map_data[8]
    a.query_sequence = map_data[9]
    a.query_qualities = pysam.qualitystring_to_array(map_data[10])

    # Set tags from tuples (tag_name, tag_value)
    for tag_name, tag_value in map_data[11:]:
        a.set_tag(tag_name, tag_value)

    return a


def map_file(
    ref_file,
    r1_file,
    r2_file,
    output_file,
    fwd_lib=True,
    max_mismatches=10,
    threads=8,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
    index_dir=None,
    index_only=False,
    batch_size=1000,
):
    """Map FASTQ reads to reference with dual-base conversion chemistry."""
    # Preflight: validate input files early to avoid spawning workers on bad paths
    if not index_only:
        if not r1_file or not os.path.exists(r1_file):
            raise FileNotFoundError(f"Input R1 file not found: {r1_file}")
        if r2_file is not None and not os.path.exists(r2_file):
            raise FileNotFoundError(f"Input R2 file not found: {r2_file}")
    if ref_file and not os.path.exists(ref_file):
        raise FileNotFoundError(f"Reference file not found: {ref_file}")
    # Determine index directory path and auto-clean temp via atexit
    if index_dir:
        os.makedirs(index_dir, exist_ok=True)
        index_base_dir = index_dir
    else:
        index_base_dir = tempfile.mkdtemp(prefix="coralsnake_")
        atexit.register(
            lambda p=index_base_dir: os.path.isdir(p)
            and shutil.rmtree(p, ignore_errors=True)
        )

    # Indexing with step-level progress (orig + MK)
    orig_fa_path = os.path.join(index_base_dir, "ref.orig.fa")
    mk_prefix_path = os.path.join(index_base_dir, "ref.mk")
    mk_index_ready = all(os.path.exists(mk_prefix_path + ext) for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"])
    orig_ready = os.path.exists(orig_fa_path)

    with Progress(
        SpinnerColumn(style="cyan"),
        TextColumn("[bold green]Index[/bold green]"),
        TextColumn("{task.description}"),
        TextColumn("| [yellow]{task.fields[conv]}[/yellow] | ORIG: [cyan]{task.fields[p1]:>5.1f}%[/cyan] | MK: [magenta]{task.fields[p2]:>5.1f}%[/magenta]"),
        transient=True,
    ) as progress:
        if index_only:
            if not ref_file:
                raise RuntimeError("--index-only requires --ref-file to be provided")
            task = progress.add_task(
                "Preparing indices...",
                total=None,
                conv="Starting",
                p1=0.0,
                p2=0.0,
            )
            def on_update(conv, p1, p2):
                progress.update(task, conv=conv, p1=p1, p2=p2)
                progress.refresh()
            _build_indices_with_progress(ref_file, index_base_dir, on_update)
            progress.update(task, description="✓ Indices ready", conv="Done", p1=100.0, p2=100.0)
            progress.refresh()
            return
        else:
            if orig_ready and mk_index_ready:
                task = progress.add_task(
                    "✓ Indices ready",
                    total=None,
                    conv="Done",
                    p1=100.0,
                    p2=100.0,
                )
                idx0_file, idx_mk_file = orig_fa_path, mk_prefix_path
            else:
                task = progress.add_task(
                    "Preparing indices...",
                    total=None,
                    conv="Starting",
                    p1=0.0,
                    p2=0.0,
                )
                def on_update(conv, p1, p2):
                    progress.update(task, conv=conv, p1=p1, p2=p2)
                    progress.refresh()
                idx0_file, idx_mk_file = _build_indices_with_progress(
                    ref_file, index_base_dir, on_update
                )
                progress.update(task, description="✓ Indices ready", conv="Done", p1=100.0, p2=100.0)
                progress.refresh()

    if index_only:
        return

    # BAM header
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": []}
    fa_for_header = (
        ref_file if ref_file else os.path.join(index_base_dir, "ref.orig.fa")
    )
    for rec in fastx_read(fa_for_header):
        header["SQ"].append({"SN": rec.name, "LN": rec.length})

    # Streaming batches in parallel (no file splitting). threads = workers
    paired = r2_file is not None
    with ExitStack() as stack:
        mode = "w" if output_file.endswith(".sam") else "wb"
        bam_out = stack.enter_context(
            pysam.AlignmentFile(output_file, mode, header=header)
        )
        progress = stack.enter_context(
            Progress(
                SpinnerColumn(style="cyan"),
                TextColumn("[bold green]Map[/bold green]"),
                TextColumn("{task.description}"),
            )
        )
        unit = "pairs" if paired else "reads"
        task = progress.add_task(
            f"[green]0[/green] / [white]0[/white] {unit} ([magenta]0.00s[/magenta])",
            total=None,
        )
        processed_reads = 0
        mapped_reads = 0

        def write_mapped(mapped):
            nonlocal processed_reads, mapped_reads
            for read_result in mapped:
                for i, item in enumerate(read_result):
                    if paired and len(item) == 3:
                        map1, map2 = item[1], item[2]
                        a1 = create_bam_record(
                            bam_out.header, map1, is_secondary=(i > 0)
                        )
                        a2 = create_bam_record(
                            bam_out.header, map2, is_secondary=(i > 0)
                        )
                        bam_out.write(a1)
                        bam_out.write(a2)
                    else:
                        map1 = item[1]
                        a1 = create_bam_record(
                            bam_out.header, map1, is_secondary=(i > 0)
                        )
                        bam_out.write(a1)
                processed_reads += 1
                if read_result:
                    mapped_reads += 1
            elapsed = format_duration(progress.tasks[task].elapsed)
            progress.update(
                task,
                description=f"[green]{mapped_reads:,}[/green] / [white]{processed_reads:,}[/white] {unit} ([magenta]{elapsed}[/magenta])",
            )

        if not paired:
            it1 = ((rec.name, rec.sequence, rec.quality) for rec in fastx_read(r1_file))
            batch = []
            with ProcessPoolExecutor(
                max_workers=max(1, threads),
                mp_context=mp.get_context("spawn"),
                initializer=_init_worker,
                initargs=(idx0_file, idx_mk_file, threads),
            ) as ex:
                futures = []
                for rec in it1:
                    batch.append(rec)
                    if len(batch) >= batch_size:
                        futures.append(
                            ex.submit(
                                _map_batch_worker,
                                list(batch),
                                False,
                                fwd_lib,
                                max_mismatches,
                                min_alignment_length,
                                min_mapping_ratio,
                            )
                        )
                        batch.clear()
                    # Drain completed futures to avoid initial stall
                    for fut in futures[:]:
                        if fut.done():
                            write_mapped(fut.result())
                            futures.remove(fut)
                if batch:
                    futures.append(
                        ex.submit(
                            _map_batch_worker,
                            list(batch),
                            False,
                            fwd_lib,
                            max_mismatches,
                            min_alignment_length,
                            min_mapping_ratio,
                        )
                    )
                for fut in as_completed(futures):
                    write_mapped(fut.result())
        else:
            it_pairs = (
                ((r1.name, r1.sequence, r1.quality), (r2.name, r2.sequence, r2.quality))
                for r1, r2 in read_paired_fastx(r1_file, r2_file)
            )
            batch = []
            with ProcessPoolExecutor(
                max_workers=max(1, threads),
                mp_context=mp.get_context("spawn"),
                initializer=_init_worker,
                initargs=(idx0_file, idx_mk_file, threads),
            ) as ex:
                futures = []
                for rec1, rec2 in it_pairs:
                    base1 = rec1[0].split()[0].rstrip("/1").rstrip("/2")
                    base2 = rec2[0].split()[0].rstrip("/1").rstrip("/2")
                    if base1 != base2:
                        raise ValueError(
                            f"r1 and r2 not in the same order: {rec1[0]} vs {rec2[0]}"
                        )
                    batch.append((rec1, rec2))
                    if len(batch) >= batch_size:
                        futures.append(
                            ex.submit(
                                _map_batch_worker,
                                list(batch),
                                True,
                                fwd_lib,
                                max_mismatches,
                                min_alignment_length,
                                min_mapping_ratio,
                            )
                        )
                        batch.clear()
                    # Drain completed futures to avoid initial stall
                    for fut in futures[:]:
                        if fut.done():
                            write_mapped(fut.result())
                            futures.remove(fut)
                if batch:
                    futures.append(
                        ex.submit(
                            _map_batch_worker,
                            list(batch),
                            True,
                            fwd_lib,
                            max_mismatches,
                            min_alignment_length,
                            min_mapping_ratio,
                        )
                    )
                for fut in as_completed(futures):
                    write_mapped(fut.result())
