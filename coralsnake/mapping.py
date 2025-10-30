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
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack

import pysam
from bwamem import BwaAligner, BwaIndexer, fastx_read, read_paired_fastx
from rich.progress import Progress, SpinnerColumn, TextColumn

from . import seqops
from .utils import convert_file_realtime, format_duration, km_conversion, mk_conversion, reverse_complement


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
    # MK (A→G, C→T) - single converted reference for both orientations
    mk_fa = os.path.join(index_base_dir, "ref.mk.fa")
    mk_prefix = os.path.splitext(mk_fa)[0]

    # Start ORIG index in parallel with conversion (they're independent)
    indexer1 = BwaIndexer(verbose=2, capture_progress=True)
    indexer2 = BwaIndexer(verbose=2, capture_progress=True)

    def build_orig():
        indexer1.build_index(orig_fa, prefix=orig_prefix)

    def build_mk():
        indexer2.build_index(mk_fa, prefix=mk_prefix)

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

        # Poll both indices (do-while pattern: always run at least once)
        while True:
            # If future is done but progress still 0, it completed too fast - set to 100%
            p1 = 100.0 if fut_orig.done() else (indexer1.progress_percent or 0.0)
            p2 = 100.0 if fut_mk.done() else (indexer2.progress_percent or 0.0)
            on_update(conv_status, p1, p2)

            if fut_orig.done() and fut_mk.done():
                break

            time.sleep(poll_interval)

        fut_orig.result()
        fut_mk.result()

    on_update("Done", 100.0, 100.0)
    return orig_fa, mk_prefix


def _map_batch_worker(
    batch,
    paired,
    max_mismatches,
    min_alignment_length,
    min_mapping_ratio,
):
    """Map a batch of reads; return (read_info, mapping_result) tuples."""
    # Use per-process cached resources initialized by _init_worker
    results = []
    if not paired:
        for name1, seq1, qua1 in batch:
            mapping_result = run_mapping_se(
                name1,
                seq1,
                qua1,
                _FORWARD_LIBRARY,
                max_mismatches,
                min_alignment_length,
                min_mapping_ratio,
            )
            # Return (read_info, mapping_result) tuple
            results.append(((name1, seq1, qua1), mapping_result))
    else:
        for (name1, seq1, qua1), (name2, seq2, qua2) in batch:
            base1 = name1.split()[0].rstrip("/1").rstrip("/2")
            base2 = name2.split()[0].rstrip("/1").rstrip("/2")
            if base1 != base2:
                raise ValueError(f"r1 and r2 not in the same order: {name1} vs {name2}")
            mapping_result = run_mapping_pe(
                name1,
                seq1,
                seq2,
                qua1,
                qua2,
                _FORWARD_LIBRARY,
                max_mismatches,
                min_alignment_length,
                min_mapping_ratio,
            )
            # Return (read_info, mapping_result) tuple
            results.append(((name1, seq1, qua1, name2, seq2, qua2), mapping_result))

    return results


# Per-process cached aligners (initialized once per worker)
_ALIGNER_ORIG = None
_ALIGNER_MK = None
_ORIENTATION_FILTER = None
_FORWARD_LIBRARY = None


def _init_worker(orig_fa, mk_index_prefix, orientation_filter, forward_library):
    global _ALIGNER_ORIG, _ALIGNER_MK, _ORIENTATION_FILTER, _FORWARD_LIBRARY
    _ORIENTATION_FILTER = orientation_filter
    _FORWARD_LIBRARY = forward_library

    # Use shorter seeds (14 instead of default 19) for better sensitivity with modified bases
    # Use moderately higher max_occ (1000 vs default 500) for repetitive elements
    # Note: The optind bug in bwamem has been fixed, so aligner creation order doesn't matter
    orig_prefix = os.path.splitext(orig_fa)[0]
    _ALIGNER_ORIG = BwaAligner(
        orig_prefix,
        min_seed_len=14,
        max_occ=1000,
        softclip_supplementary=True,
        mark_secondary=True,
        clip_penalties=(6, 6),
        unpaired_penalty=24,
        min_score=20,
        insert_model=(80, 60, 450),
    )

    _ALIGNER_MK = BwaAligner(
        mk_index_prefix,
        min_seed_len=14,
        max_occ=1000,
        softclip_supplementary=True,
        mark_secondary=True,
        clip_penalties=(6, 6),
        unpaired_penalty=24,
        min_score=20,
        insert_model=(80, 60, 450),
    )


def find_properly_paired_hits(hits, fwd=True):
    """Find read1/read2 hit pairs on the same contig, within 1 kb.
    
    Note: Since we pre-convert Read2 with RC before alignment, both reads may map
    to the same strand. The original opposite-strand check has been removed.
    """
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
                    # Removed opposite-strand check since we pre-RC Read2
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
    forward_library=True,
    max_mismatches=10,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
):
    """Map one single-end read and return scored alignments."""
    mapped = []
    # Filter orientations if specified
    orientations = [1, 2] if _ORIENTATION_FILTER is None else [_ORIENTATION_FILTER]
    for orientation in orientations:
        # Build converted read
        if orientation == 1:
            seq1_conv = mk_conversion(seq1) if forward_library else km_conversion(seq1)
        else:
            seq1_conv = km_conversion(seq1) if forward_library else mk_conversion(seq1)

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
                s = seqops.reverse_complement(seq1)
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
    forward_library=True,
    max_mismatches=10,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
):
    """Map one paired-end read and return scored alignment pairs."""
    mapped = []
    # Filter orientations if specified
    orientations = [1, 2] if _ORIENTATION_FILTER is None else [_ORIENTATION_FILTER]
    for orientation in orientations:
        # Build converted reads
        # CRITICAL: Both reads MUST use the same chemistry (MK) for mate rescue to work!
        # Since PE reads sequence opposite strands, Read2 needs RC before MK conversion.
        if orientation == 1:
            if forward_library:
                # Orientation 1, forward library
                # Read1: forward strand → MK conversion
                # Read2: reverse strand → RC + MK conversion
                seq1_conv = mk_conversion(seq1)
                seq2_conv = mk_conversion(reverse_complement(seq2))
            else:
                # Orientation 1, reverse library
                seq1_conv = mk_conversion(reverse_complement(seq1))
                seq2_conv = mk_conversion(seq2)
        else:
            if forward_library:
                # Orientation 2, forward library (opposite of orientation 1)
                seq1_conv = mk_conversion(reverse_complement(seq1))
                seq2_conv = mk_conversion(seq2)
            else:
                # Orientation 2, reverse library
                seq1_conv = mk_conversion(seq1)
                seq2_conv = mk_conversion(reverse_complement(seq2))

        # Align reads using BWA's paired-end mode (with mate rescue)
        pe_alignments = tuple(_ALIGNER_MK.align(seq1_conv, seq2_conv))
        
        # Extract hits from PE alignments and mark read numbers
        # pe_alignments is a tuple of PairedAlignment(read1, read2, is_proper_pair, insert_size)
        # Deduplicate hits by (ctg, r_st, r_en, strand) since PE mode can return
        # the same hit multiple times in different pair combinations
        seen_hits1 = {}
        seen_hits2 = {}
        
        for paired_aln in pe_alignments:
            if paired_aln.read1:
                hit_key = (paired_aln.read1.ctg, paired_aln.read1.r_st, paired_aln.read1.r_en, paired_aln.read1.strand)
                if hit_key not in seen_hits1:
                    paired_aln.read1.read_num = 1
                    seen_hits1[hit_key] = paired_aln.read1
            if paired_aln.read2:
                hit_key = (paired_aln.read2.ctg, paired_aln.read2.r_st, paired_aln.read2.r_en, paired_aln.read2.strand)
                if hit_key not in seen_hits2:
                    paired_aln.read2.read_num = 2
                    seen_hits2[hit_key] = paired_aln.read2
        
        hits1 = list(seen_hits1.values())
        hits2 = list(seen_hits2.values())
        combined_hits = tuple(hits1) + tuple(hits2)
        
        filtered_hits = filter_hits(
            combined_hits,
            seq1,
            seq2,
            min_alignment_length,
            min_mapping_ratio,
        )
        
        # Try to find properly paired hits (same contig)
        paired_hits = list(find_properly_paired_hits(filtered_hits, fwd=True))
        
        # If no same-contig pairs found, try cross-contig rescue
        if not paired_hits:
            # Separate hits by read number
            read1_hits = [h for h in filtered_hits if hasattr(h, 'read_num') and h.read_num == 1]
            read2_hits = [h for h in filtered_hits if hasattr(h, 'read_num') and h.read_num == 2]
            
            if read1_hits and read2_hits:
                # Try rescue: use best hits from each read as anchors
                # For each Read1 hit, check if any Read2 hit could pair with it (and vice versa)
                # This allows cross-contig pairing which mate rescue might have found
                for h1 in read1_hits[:3]:  # Try top 3 Read1 hits as anchors
                    for h2 in read2_hits[:3]:  # Try top 3 Read2 hits as anchors
                        # If reads are on different contigs, consider them as potential rescued pairs
                        # This handles cases where mem_matesw found a hit on a different contig
                        if h1.ctg != h2.ctg:
                            # Cross-contig pair - likely from mate rescue
                            paired_hits.append((h1, h2))
        
        for hit1, hit2 in paired_hits:
            tlen = max(hit1.r_en, hit2.r_en) - min(hit1.r_st, hit2.r_st)
            # Fetch original reference for scoring
            ref1 = _ALIGNER_ORIG.seq(hit1.ctg, hit1.r_st, hit1.r_en)
            ref2 = _ALIGNER_ORIG.seq(hit2.ctg, hit2.r_st, hit2.r_en)
            
            read1_reverse = hit1.strand == -1
            read2_reverse = hit2.strand == -1
            
            # For scoring, we need to account for the fact that Read2 was RC'd before conversion
            # Since we did: align(MK(seq1), MK(RC(seq2))) to ref_mk
            # And hits map to forward strand, we need to use RC(seq2) for scoring
            if read1_reverse:
                s1 = seqops.reverse_complement(seq1)
                q1 = qua1[::-1]
            else:
                s1 = seq1
                q1 = qua1
                
            # CRITICAL: Read2 was RC'd before alignment for orientation 1
            # So we must use RC(seq2) for scoring, regardless of hit2.strand
            if orientation == 1:
                if forward_library:
                    # seq2 was RC'd before conversion, so use RC for scoring
                    s2 = seqops.reverse_complement(seq2)
                    q2 = qua2[::-1]
                else:
                    s2 = seq2
                    q2 = qua2
            else:  # orientation == 2
                if forward_library:
                    s2 = seq2
                    q2 = qua2
                else:
                    s2 = seqops.reverse_complement(seq2)
                    q2 = qua2[::-1]

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
        
        # If no paired hits found, keep single-read mappings (mate unmapped)
        if not paired_hits:
            # Separate hits by read number
            read1_hits = [h for h in filtered_hits if hasattr(h, 'read_num') and h.read_num == 1]
            read2_hits = [h for h in filtered_hits if hasattr(h, 'read_num') and h.read_num == 2]
            
            # Process read1 single mappings
            for hit in read1_hits:
                ref = _ALIGNER_ORIG.seq(hit.ctg, hit.r_st, hit.r_en)
                read_reverse = hit.strand == -1
                if read_reverse:
                    # Flag 89 = paired, mapped, mate unmapped, read reverse, first in pair
                    flag = 89
                    s = seqops.reverse_complement(seq1)
                    q = qua1[::-1]
                else:
                    # Flag 73 = paired, mapped, mate unmapped, first in pair  
                    flag = 73
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
                    "*",  # Mate unmapped
                    0,
                    0,
                    s,
                    q,
                ] + tags
                mapped.append([score, map1])
            
            # Process read2 single mappings
            for hit in read2_hits:
                ref = _ALIGNER_ORIG.seq(hit.ctg, hit.r_st, hit.r_en)
                read_reverse = hit.strand == -1
                if read_reverse:
                    # Flag 153 = paired, mapped, mate unmapped, read reverse, second in pair
                    flag = 153
                    s = seqops.reverse_complement(seq2)
                    q = qua2[::-1]
                else:
                    # Flag 137 = paired, mapped, mate unmapped, second in pair
                    flag = 137
                    s = seq2
                    q = qua2
                
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
                map2 = [
                    name,
                    flag,
                    hit.ctg,
                    hit.r_st + 1,
                    mapq,
                    cigar_str,
                    "*",  # Mate unmapped
                    0,
                    0,
                    s,
                    q,
                ] + tags
                mapped.append([score, map2])

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


def create_unmapped_record(header, name, seq, qual, flag):
    """Create an unmapped pysam.AlignedSegment."""
    a = pysam.AlignedSegment(header=header)
    a.query_name = name.split()[0]
    a.flag = flag
    a.reference_id = -1
    a.reference_start = -1
    a.mapping_quality = 0
    a.cigarstring = None
    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = 0
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(qual)
    return a


def map_file(
    ref_file,
    r1_file,
    r2_file,
    output_file,
    forward_library=True,
    max_mismatches=10,
    threads=8,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
    index_dir=None,
    index_only=False,
    batch_size=1000,
    orientation_filter=None,
):
    """Map FASTQ reads to reference with dual-base conversion chemistry.

    Args:
        forward_library: True for forward library, False for reverse library.
        orientation_filter: If specified, only map to this orientation (1 or 2).
                          None means map to both orientations (default).
    """
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
    mk_index_ready = all(
        os.path.exists(mk_prefix_path + ext)
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]
    )
    orig_ready = os.path.exists(orig_fa_path)

    with Progress(
        SpinnerColumn(style="cyan"),
        TextColumn("[bold green]Index[/bold green]"),
        TextColumn("{task.description}"),
        TextColumn(
            "| [yellow]{task.fields[conv]}[/yellow] | ORIG: [cyan]{task.fields[p1]:>5.1f}%[/cyan] | MK: [cyan]{task.fields[p2]:>5.1f}%[/cyan] ([magenta]{task.fields[elapsed]}[/magenta])"
        ),
        transient=False,
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
                elapsed="0.00s",
            )
            start_time = time.time()

            def on_update(conv, p1, p2):
                elapsed = time.time() - start_time
                elapsed_str = format_duration(elapsed)
                progress.update(task, conv=conv, p1=p1, p2=p2, elapsed=elapsed_str)
                progress.refresh()
                time.sleep(0.001)  # Force thread context switch for Rich display

            _build_indices_with_progress(ref_file, index_base_dir, on_update)
            elapsed = time.time() - start_time
            progress.update(
                task,
                description="✓ Indices ready",
                conv="Done",
                p1=100.0,
                p2=100.0,
                elapsed=format_duration(elapsed),
            )
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
                    elapsed="0.00s",
                )
                idx0_file, idx_mk_file = orig_fa_path, mk_prefix_path
            else:
                task = progress.add_task(
                    "Preparing indices...",
                    total=None,
                    conv="Starting",
                    p1=0.0,
                    p2=0.0,
                    elapsed="0.00s",
                )
                start_time = time.time()

                def on_update(conv, p1, p2):
                    elapsed = time.time() - start_time
                    elapsed_str = format_duration(elapsed)
                    progress.update(task, conv=conv, p1=p1, p2=p2, elapsed=elapsed_str)
                    progress.refresh()
                    time.sleep(0.001)  # Force thread context switch for Rich display

                idx0_file, idx_mk_file = _build_indices_with_progress(
                    ref_file, index_base_dir, on_update
                )
                elapsed = time.time() - start_time
                progress.update(
                    task,
                    description="✓ Indices ready",
                    conv="Done",
                    p1=100.0,
                    p2=100.0,
                    elapsed=format_duration(elapsed),
                )
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

        def write_mapped(batch_results):
            nonlocal processed_reads, mapped_reads
            for read_info, mapping_result in batch_results:
                if mapping_result:
                    # Write mapped reads
                    for i, item in enumerate(mapping_result):
                        if paired and len(item) == 3:
                            # Properly paired reads
                            map1, map2 = item[1], item[2]
                            a1 = create_bam_record(
                                bam_out.header, map1, is_secondary=(i > 0)
                            )
                            a2 = create_bam_record(
                                bam_out.header, map2, is_secondary=(i > 0)
                            )
                            bam_out.write(a1)
                            bam_out.write(a2)
                        elif paired and len(item) == 2:
                            # Single-read mapping from PE mode (mate unmapped)
                            map1 = item[1]
                            a1 = create_bam_record(
                                bam_out.header, map1, is_secondary=(i > 0)
                            )
                            bam_out.write(a1)
                            # Also write the unmapped mate
                            name1, seq1, qua1, name2, seq2, qua2 = read_info
                            # Determine which read is mapped based on flag
                            if map1[1] & 0x40:  # First in pair bit (0x40)
                                # Read1 is mapped, write unmapped read2
                                # Base flag: 128 (read2) + 4 (unmapped) + 1 (paired) = 133
                                # Add 32 (mate reverse) if read1 is reverse
                                flag2 = 133
                                if map1[1] & 0x10:  # Mapped mate is reverse
                                    flag2 |= 0x20  # Set mate reverse bit
                                a2 = create_unmapped_record(
                                    bam_out.header, name2, seq2, qua2, flag2
                                )
                                a2.next_reference_name = map1[2]  # Mate's contig
                                a2.next_reference_start = map1[3] - 1  # Mate's position
                                bam_out.write(a2)
                            else:  # Second in pair (0x80)
                                # Read2 is mapped, write unmapped read1
                                # Base flag: 64 (read1) + 4 (unmapped) + 1 (paired) = 69
                                # Add 32 (mate reverse) if read2 is reverse
                                flag1 = 69
                                if map1[1] & 0x10:  # Mapped mate is reverse
                                    flag1 |= 0x20  # Set mate reverse bit
                                a1_unmapped = create_unmapped_record(
                                    bam_out.header, name1, seq1, qua1, flag1
                                )
                                a1_unmapped.next_reference_name = map1[2]  # Mate's contig
                                a1_unmapped.next_reference_start = map1[3] - 1  # Mate's position
                                bam_out.write(a1_unmapped)
                        else:
                            # Single-end mapping
                            map1 = item[1]
                            a1 = create_bam_record(
                                bam_out.header, map1, is_secondary=(i > 0)
                            )
                            bam_out.write(a1)
                    mapped_reads += 1
                else:
                    # Write unmapped reads
                    if paired:
                        name1, seq1, qua1, name2, seq2, qua2 = read_info
                        # Flag 77 = paired, unmapped, mate unmapped, first in pair
                        a1 = create_unmapped_record(
                            bam_out.header, name1, seq1, qua1, 77
                        )
                        # Flag 141 = paired, unmapped, mate unmapped, second in pair
                        a2 = create_unmapped_record(
                            bam_out.header, name2, seq2, qua2, 141
                        )
                        bam_out.write(a1)
                        bam_out.write(a2)
                    else:
                        name1, seq1, qua1 = read_info
                        # Flag 4 = unmapped
                        a1 = create_unmapped_record(
                            bam_out.header, name1, seq1, qua1, 4
                        )
                        bam_out.write(a1)
                processed_reads += 1
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
                initargs=(idx0_file, idx_mk_file, orientation_filter, forward_library),
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
                initargs=(idx0_file, idx_mk_file, orientation_filter, forward_library),
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
                            max_mismatches,
                            min_alignment_length,
                            min_mapping_ratio,
                        )
                    )
                for fut in as_completed(futures):
                    write_mapped(fut.result())
