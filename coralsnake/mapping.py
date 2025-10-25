#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-08 20:32


import os
import random
import tempfile
import shutil

import mappy as mp
import pysam
from rich.progress import Progress

from .utils import format_duration, mk_conversion, km_conversion, convert_file_realtime
from . import seqops


def split_fastq_files(r1_file, r2_file, num_chunks, output_dir):
    """
    Split FASTQ file(s) into chunks for parallel processing using fast C implementation.
    
    Supports both .fq and .fq.gz files. Automatically estimates reads per chunk based on file size.
    
    Args:
        r1_file: Path to read 1 FASTQ file (.fq or .fq.gz)
        r2_file: Path to read 2 FASTQ file (None for single-end, .fq or .fq.gz)
        num_chunks: Number of chunks to split into
        output_dir: Directory to write chunk files
        
    Returns:
        List of tuples: [(r1_chunk1, r2_chunk1), (r1_chunk2, r2_chunk2), ...]
        For single-end: [(r1_chunk1, None), (r1_chunk2, None), ...]
    """
    # Split R1 file using fast C implementation
    r1_chunks = seqops.fast_split_fastq(r1_file, output_dir, num_chunks, False)
    
    if r2_file:
        # Split R2 file for paired-end
        r2_chunks = seqops.fast_split_fastq(r2_file, output_dir, num_chunks, True)
        # Pair up R1 and R2 chunks
        chunk_files = list(zip(r1_chunks, r2_chunks))
    else:
        # Single-end: pair with None
        chunk_files = [(r1_chunk, None) for r1_chunk in r1_chunks]
    
    return chunk_files


def process_chunk(
    chunk_r1,
    chunk_r2,
    chunk_output,
    idx0_file,
    idx_mk_file,
    ref_file,
    fwd_lib,
    max_mismatches,
    threads,
    min_alignment_length,
    min_mapping_ratio,
    chunk_id,
):
    """
    Process a single chunk of FASTQ files.
    
    This function is designed to be called for each chunk.
    """
    # Load indices
    idx0 = mp.Aligner(fn_idx_in=idx0_file, preset="sr", n_threads=threads)
    idx_mk = mp.Aligner(fn_idx_in=idx_mk_file, preset="sr", n_threads=threads)
    
    # Create BAM header
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": []}
    for name, seq, *_ in mp.fastx_read(ref_file):
        header["SQ"].append({"SN": name, "LN": len(seq)})
    
    # Process reads and write to BAM
    with pysam.AlignmentFile(chunk_output, "wb", header=header) as bam_out:
        if chunk_r2 is None:
            # Single-end
            for name1, seq1, qua1 in mp.fastx_read(chunk_r1):
                mapped = run_mapping(
                    name1,
                    seq1,
                    None,
                    qua1,
                    None,
                    idx0,
                    idx_mk,
                    fwd_lib,
                    max_mismatches,
                    min_alignment_length,
                    min_mapping_ratio,
                )
                
                for i, item in enumerate(mapped):
                    map1 = item[1]
                    a1 = create_bam_record(bam_out.header, map1, is_secondary=(i > 0))
                    bam_out.write(a1)
        else:
            # Paired-end
            for (name1, seq1, qua1), (name2, seq2, qua2) in zip(
                mp.fastx_read(chunk_r1), mp.fastx_read(chunk_r2)
            ):
                if name1 != name2:
                    raise ValueError("r1 and r2 not in the same order")
                
                mapped = run_mapping(
                    name1,
                    seq1,
                    seq2,
                    qua1,
                    qua2,
                    idx0,
                    idx_mk,
                    fwd_lib,
                    max_mismatches,
                    min_alignment_length,
                    min_mapping_ratio,
                )
                
                for i, item in enumerate(mapped):
                    map1, map2 = item[1], item[2]
                    a1 = create_bam_record(bam_out.header, map1, is_secondary=(i > 0))
                    a2 = create_bam_record(bam_out.header, map2, is_secondary=(i > 0))
                    bam_out.write(a1)
                    bam_out.write(a2)
    
    return chunk_output


def find_properly_paired_hits(hits, fwd=True):
    """
    # find properly paired hits
    Given a list of hits with read_num, ref_name, ref_start, ref_end, strand. For example:
        1 rRNA-Hsa-nucleus_locus 8204 8237 1
        1 rRNA-Hsa-nucleus_locus 21555 21588 1
        2 rRNA-Hsa-nucleus_locus 21561 21596 -1
        2 rRNA-Hsa-nucleus_locus 8210 8245 -1

    This function will find all properly paired hits. A hit is properly paired if it is properly paired.
    We defined properly paired hits as:
    1. The hit paris are from two reads (read_num=1 and read_num=2)
    2. Both hits mapped to the same reference sequence
    3. read 1 and reads 2 of one pair are in different strand (if reads 1 is 1, reads 2 is -1)
    4. If fwd is True, strand = 1 is forward, strand = -1 is reverse. Otherwise, strand = 1 is reverse, strand = -1 is forward
    5. The mapping start of read 1 is smaller than the mapping end of read 2, if fwd is True. Otherwise, the mapping end of read 1 is larger than the mapping start of read 2.
    6. The distance between the mapping start and end of read 1 and read 2 is less than 1000
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
                    if hit1.strand + hit2.strand == 0:
                        if fwd:
                            if hit1.r_st < hit2.r_en and hit2.r_en - hit1.r_st < 1000:
                                parsed_hits.append((hit1, hit2))
                        else:
                            if hit1.r_en > hit2.r_st and hit1.r_en - hit2.r_st < 1000:
                                parsed_hits.append((hit1, hit2))

    return parsed_hits


def cal_md_and_tag(cigar, seq, ref, fwd):
    """
    Calculate MD tag and custom tags for dual-base conversion chemistry.

    This function generates the MD tag (mismatch/deletion string) and calculates
    custom tags for tracking conversion statistics in dual-base conversion sequencing.

    Dual-base conversion chemistry (not standard bisulfite):
    - MK conversion: C->T AND A->G simultaneously
    - KM conversion: G->A AND T->C simultaneously

    Args:
        cigar: List of CIGAR operations [(length, operation), ...]
        seq: Query sequence (read)
        ref: Reference sequence
        fwd: If True, expect MK conversions (C->T, A->G)
             If False, expect KM conversions (G->A, T->C)

    Returns:
        Tuple of (md_tag, yf, zf, yc, zc, ns, nc) where:
        - md_tag: MD tag string for SAM/BAM format
        - yf: Number of A->G conversions (MK) or T->C conversions (KM)
        - zf: Number of A->A matches (MK) or T->T matches (KM)
        - yc: Number of C->T conversions (MK) or G->A conversions (KM)
        - zc: Number of C->C matches (MK) or G->G matches (KM)
        - ns: Number of non-conversion mismatches (sequencing errors)
        - nc: Number of clipped bases and indels

    CIGAR operations:
        M (0): Alignment match/mismatch
    I (1): Insertion to the reference
    D (2): Deletion from the reference
    S (4): Soft clipping
    """
    # Use optimized C implementation
    return seqops.fast_cal_md_and_tag(cigar, seq, ref, fwd)


def calculate_directional_score(cigar, seq, ref, is_orientation1):
    """
    Calculate a directional conversion-aware alignment score for dual-base conversion chemistry.

    IMPORTANT: The conversion is NOT 100% efficient!

    For orientation 1 (MK conversion applied to read):
        - GOOD: Perfect matches (any base matches)
        - GOOD: Expected conversions C->T, A->G (chemistry worked)
        - GOOD: Unconverted C->C, A->A (chemistry didn't work, but that's OK)
        - BAD: Wrong conversions T->C, G->A (wrong chemistry applied)
        - BAD: Other mismatches (sequencing errors)

    For orientation 2 (KM conversion applied to read):
        - GOOD: Perfect matches (any base matches)
        - GOOD: Expected conversions G->A, T->C (chemistry worked)
        - GOOD: Unconverted G->G, T->T (chemistry didn't work, but that's OK)
        - BAD: Wrong conversions A->G, C->T (wrong chemistry applied)
        - BAD: Other mismatches (sequencing errors)

    The key insight: when comparing unconverted read vs unconverted reference:
        - Perfect match = GOOD
        - Expected conversion = GOOD (chemistry worked)
        - Note: Unconverted bases will show as perfect matches!
        - Wrong conversion = BAD (indicates wrong orientation)
        - Other mismatch = BAD (sequencing error)

    Returns:
        score: alignment score (matches + expected_conversions - wrong_conversions - other_mismatches - indels)
        wrong_conversions: number of wrong-direction conversions
        total_bad_mismatches: wrong_conversions + other_mismatches
    """
    # Use optimized C implementation
    return seqops.fast_calculate_directional_score(cigar, seq, ref, is_orientation1)


def filter_hits(hits, seq1, seq2, min_alignment_length=20, min_mapping_ratio=0.5):
    """
    Filter hits by the following rules:
    1. The mapping quality of the hit is greater than 0
    2. The alignment length of the hit is greater than min_alignment_length
    3. The mapping length of the hit is larger than min_mapping_ratio of the query length

    Args:
        hits: List of alignment hits
        seq1: Read 1 sequence
        seq2: Read 2 sequence (or None for single-end)
        min_alignment_length: Minimum alignment length (default: 20)
        min_mapping_ratio: Minimum mapping length ratio (default: 0.5)
    """
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


def run_mapping(
    name,
    seq1,
    seq2,
    qua1,
    qua2,
    idx0,
    idx_mk,
    fwd_lib=True,
    max_mismatches=10,
    min_alignment_length=20,
    min_mapping_ratio=0.5,
):
    """
    Map reads using dual-base conversion chemistry with directional filtering.

    This function implements a directional mapping strategy for dual-base conversion
    sequencing (similar to bisulfite-seq but converts TWO bases simultaneously).
    It tries both conversion orientations and filters based on conversion patterns.

    Conversion chemistry:
        - MK conversion: C->T AND A->G simultaneously
        - KM conversion: G->A AND T->C simultaneously (reverse complement)

    Mapping strategy:
        - Uses only MK converted reference
        - Tries two orientations:
            * Orientation 1: read1 MK, read2 KM (for forward library)
            * Orientation 2: read1 KM, read2 MK (for forward library)
        - Filters alignments based on wrong-direction conversions

    Args:
        name: Read name/identifier
        seq1: Read 1 sequence
        seq2: Read 2 sequence (None for single-end)
        qua1: Read 1 quality string
        qua2: Read 2 quality string (None for single-end)
        idx0: Original (unconverted) reference index
        idx_mk: MK converted reference index
        fwd_lib: True for forward library, False for reverse library
        max_mismatches: Maximum allowed bad mismatches (wrong conversions + errors)

    Returns:
        List of [score, map1, map2] for paired-end or [score, map1] for single-end
        where map1/map2 are lists of SAM fields plus custom tags
    """
    idx = idx_mk

    mapped = []
    # Try both orientations
    for orientation in [1, 2]:
        if orientation == 1:
            # Orientation 1: MK conversion
            if fwd_lib:
                seq1_conv = mk_conversion(seq1)
                seq2_conv = km_conversion(seq2) if seq2 else None
            else:
                # For reverse library, read1 behaves like read2, read2 behaves like read1
                seq1_conv = km_conversion(seq1)
                seq2_conv = mk_conversion(seq2) if seq2 else None
        else:
            # Orientation 2: KM conversion
            if fwd_lib:
                seq1_conv = km_conversion(seq1)
                seq2_conv = mk_conversion(seq2) if seq2 else None
            else:
                # For reverse library, read1 behaves like read2, read2 behaves like read1
                seq1_conv = mk_conversion(seq1)
                seq2_conv = km_conversion(seq2) if seq2 else None

        if seq2:
            # Paired-end mapping
            for hit1, hit2 in find_properly_paired_hits(
                filter_hits(
                    idx.map(seq1_conv, seq2=seq2_conv, cs=False, MD=False),
                    seq1,
                    seq2,
                    min_alignment_length,
                    min_mapping_ratio,
                ),
                fwd=True,  # Always use forward for MK reference
            ):
                # Filter by strand constraints
                # Note: We don't filter by strand here - the directional score will handle it
                # The strand just tells us which strand of the reference the read mapped to
                tlen = max(hit1.r_en, hit2.r_en) - min(hit1.r_st, hit2.r_st)
                # Get reference sequences from original (unconverted) reference
                ref1 = idx0.seq(hit1.ctg, hit1.r_st, hit1.r_en)
                ref2 = idx0.seq(hit2.ctg, hit2.r_st, hit2.r_en)

                # Determine sequences and flags based on hit strand
                # hit.strand tells us the strand relative to the reference
                read1_reverse = hit1.strand == -1
                read2_reverse = hit2.strand == -1

                if read1_reverse:
                    s1 = mp.revcomp(seq1)
                    q1 = qua1[::-1]
                else:
                    s1 = seq1
                    q1 = qua1

                if read2_reverse:
                    s2 = mp.revcomp(seq2)
                    q2 = qua2[::-1]
                else:
                    s2 = seq2
                    q2 = qua2

                # Set flags based on strand
                if read1_reverse and not read2_reverse:
                    flag1, flag2 = 83, 163  # read1 reverse, read2 forward
                elif not read1_reverse and read2_reverse:
                    flag1, flag2 = 99, 147  # read1 forward, read2 reverse
                else:
                    # Both on same strand - shouldn't happen for proper pairs
                    flag1, flag2 = 67, 131  # both forward

                # fix soft clip
                # https://github.com/lh3/minimap2/issues/356
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

                # Calculate directional score to filter wrong-direction conversions
                is_orientation1 = orientation == 1
                score1, wrong_conv1, bad_mm1 = calculate_directional_score(
                    cigar1, s1, ref1, is_orientation1
                )
                score2, wrong_conv2, bad_mm2 = calculate_directional_score(
                    cigar2, s2, ref2, is_orientation1
                )

                # Filter: skip if too many bad mismatches (wrong conversions + other errors)
                # Allow some sequencing errors and incomplete conversion
                total_bad = bad_mm1 + bad_mm2

                if total_bad > max_mismatches:
                    continue

                # Calculate MD tag and custom tags for both reads
                md1, yf1, zf1, yc1, zc1, ns1, nc1 = cal_md_and_tag(
                    cigar1, s1, ref1, is_orientation1
                )
                md2, yf2, zf2, yc2, zc2, ns2, nc2 = cal_md_and_tag(
                    cigar2, s2, ref2, is_orientation1
                )

                combined_score = score1 + score2

                # Calculate MAPQ based on alignment score (simple formula: min(60, score))
                # For paired-end, use the minimum of the two scores
                mapq = min(60, min(score1, score2))

                # Create common tags for both reads
                common_tags = [
                    ("ST", orientation),
                ]

                # Create map data for read1
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
                map1 = create_map_data(
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
                    tags1,
                )

                # Create map data for read2
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
                map2 = create_map_data(
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
                    tags2,
                )

                mapped.append([combined_score, map1, map2])
        else:
            # Single-end mapping
            for hit in filter_hits(
                idx.map(seq1_conv, cs=False, MD=False),
                seq1,
                None,
                min_alignment_length,
                min_mapping_ratio,
            ):
                # Note: We don't filter by strand - the directional score will handle it
                # Get reference sequence
                ref = idx0.seq(hit.ctg, hit.r_st, hit.r_en)

                # Determine sequence and flag based on hit strand
                read_reverse = hit.strand == -1

                if read_reverse:
                    flag = 16  # reverse strand
                    s = mp.revcomp(seq1)
                    q = qua1[::-1]
                else:
                    flag = 0  # forward strand
                    s = seq1
                    q = qua1

                # Fix soft clip
                cigar_str = hit.cigar_str
                cigar = hit.cigar
                if hit.q_st > 0:
                    cigar_str = f"{hit.q_st}S" + cigar_str
                    cigar = [[hit.q_st, 4]] + cigar
                if hit.q_en < len(s):
                    cigar_str = cigar_str + f"{len(s) - hit.q_en}S"
                    cigar = cigar + [[len(s) - hit.q_en, 4]]

                # Calculate directional score to filter wrong-direction conversions
                is_orientation1 = orientation == 1
                score, wrong_conv, bad_mm = calculate_directional_score(
                    cigar, s, ref, is_orientation1
                )

                # Filter: skip if too many bad mismatches (wrong conversions + other errors)
                # Allow some sequencing errors and incomplete conversion
                # Use half the threshold for single-end reads
                if bad_mm > max_mismatches // 2:
                    continue

                # Calculate MD tag and custom tags
                md, yf, zf, yc, zc, ns, nc = cal_md_and_tag(
                    cigar, s, ref, is_orientation1
                )

                # Calculate MAPQ based on alignment score (simple formula: min(60, score))
                mapq = min(60, score)

                # Create tags
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
                map1 = create_map_data(
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
                    tags,
                )
                mapped.append([score, map1])

    random.shuffle(mapped)
    mapped = sorted(mapped, key=lambda x: x[0], reverse=True)
    return mapped


def create_map_data(
    name, flag, ctg, pos, mapq, cigar_str, mate_ctg, mate_pos, tlen, seq, qual, tags
):
    """
    Create a map data list for BAM record.

    Args:
        name: Read name
        flag: SAM flag
        ctg: Reference contig
        pos: Reference position (1-based)
        mapq: Mapping quality
        cigar_str: CIGAR string
        mate_ctg: Mate reference contig
        mate_pos: Mate reference position (1-based)
        tlen: Template length
        seq: Query sequence
        qual: Quality string
        tags: List of (tag_name, tag_value) tuples

    Returns:
        List containing SAM fields and tags
    """
    return [
        name,
        flag,
        ctg,
        pos,
        mapq,
        cigar_str,
        mate_ctg,
        mate_pos,
        tlen,
        seq,
        qual,
    ] + tags


def create_bam_record(header, map_data, is_secondary):
    """
    Create a BAM record from mapping data.

    Args:
        header: BAM header
        map_data: List containing SAM fields and tags [name, flag, ref, pos, mapq, cigar, ...]
        is_secondary: Whether this is a secondary alignment

    Returns:
        pysam.AlignedSegment object
    """
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
    num_chunks=1,
):
    """
    Map FASTQ reads to reference genome using dual-base conversion chemistry.

    This is the main entry point for the mapping pipeline. It handles reference
    conversion, index loading, read processing, and BAM file generation.

    Args:
        ref_file: Path to reference FASTA file
        r1_file: Path to read 1 FASTQ file (can be gzipped)
        r2_file: Path to read 2 FASTQ file (None for single-end, can be gzipped)
        output_file: Path to output BAM file
        fwd_lib: True for forward library, False for reverse library
        max_mismatches: Maximum allowed bad mismatches for filtering
        threads: Number of threads for minimap2 indexing (default: 8)
                 Note: Currently only used for index loading, not for parallel mapping.
                 Parallel read mapping is planned for future versions.
        min_alignment_length: Minimum alignment length for filtering (default: 20)
        min_mapping_ratio: Minimum mapping length ratio for filtering (default: 0.5)
        index_dir: Directory to store/load index files. If None, uses temporary directory (default: None)
        index_only: If True, only build indices without mapping (default: False)
        num_chunks: Number of chunks to split input files for processing (default: 1)
                    If > 1, splits FASTQ files and processes chunks sequentially, then merges BAM files

    Output BAM tags:
        - MD:Z: Mismatch/deletion string (standard SAM tag)
        - ST:i: Orientation (1=MK conversion: C->T, A->G; 2=KM conversion: G->A, T->C)
        - AS:i: Alignment score (directional conversion-aware, higher is better)
        - Yf:i: Number of forward conversions (A->G for MK, T->C for KM)
        - Zf:i: Number of forward matches (A->A for MK, T->T for KM)
        - Yc:i: Number of C conversions (C->T for MK, G->A for KM)
        - Zc:i: Number of C matches (C->C for MK, G->G for KM)
        - NS:i: Number of non-conversion mismatches (sequencing errors)
        - NC:i: Number of clipped bases and indels

    Note: Alignments are sorted by AS score (best first) for primary/secondary assignment.
    """
    # Determine index directory and cleanup behavior
    if index_dir:
        # User-specified directory: create if doesn't exist, don't clean up
        os.makedirs(index_dir, exist_ok=True)
        use_temp_dir = False
        index_base_dir = index_dir
    else:
        # Temporary directory: will be cleaned up
        use_temp_dir = True
        index_base_dir = tempfile.mkdtemp(prefix="coralsnake_")
    
    mk_file = os.path.join(index_base_dir, "ref.mk.fa")
    idx0_file = os.path.join(index_base_dir, "ref.orig.mmi")
    idx_mk_file = os.path.join(index_base_dir, "ref.mk.mmi")
    try:
        # Check if indices already exist
        indices_exist = os.path.exists(idx0_file) and os.path.exists(idx_mk_file)
        
        if indices_exist and not index_only:
            # Load existing indices
            with Progress(transient=True) as progress:
                task = progress.add_task("Loading existing indices...", total=None)
                idx0 = mp.Aligner(fn_idx_in=idx0_file, preset="sr", n_threads=threads)
                idx_mk = mp.Aligner(fn_idx_in=idx_mk_file, preset="sr", n_threads=threads)
                progress.update(task, description="✓ Indices loaded from disk")
        else:
            # Build new indices
            with Progress(transient=True) as progress:
                task = progress.add_task("Converting reference...", total=None)
                # Convert reference to MK (A->G, C->T) using fast C implementation
                convert_file_realtime(ref_file, mk_file, "AC", "GT")
                progress.update(task, description="✓ Reference converted")
                # Load original reference and save index
                progress.update(task, description="Indexing original reference...")
                idx0 = mp.Aligner(
                    fn_idx_in=ref_file,
                    preset="sr",
                    n_threads=threads,
                    k=10,
                    w=10,
                    min_cnt=0,
                    min_chain_score=0,
                    best_n=50,
                    fn_idx_out=idx0_file,
                )
                progress.update(task, description="✓ Original reference indexed")
                # Load MK converted reference and save index
                progress.update(task, description="Indexing converted reference...")
                idx_mk = mp.Aligner(
                    fn_idx_in=mk_file,
                    preset="sr",
                    n_threads=threads,
                    k=10,
                    w=10,
                    min_cnt=0,
                    min_chain_score=0,
                    best_n=50,
                    fn_idx_out=idx_mk_file,
                )
                progress.update(task, description="✓ Converted reference indexed")

                # Remove MK FASTA file after indexing (index file is sufficient)
                if os.path.exists(mk_file):
                    os.remove(mk_file)
                progress.update(task, description="✓ Index files created")
        
        # If index-only mode, we're done
        if index_only:
            return

        # Handle chunked processing if num_chunks > 1
        if num_chunks > 1:
            # Create temporary directory for chunks
            chunks_dir = tempfile.mkdtemp(prefix="coralsnake_chunks_")
            try:
                # Split FASTQ files into chunks
                print(f"\n📦 Splitting input files into {num_chunks} chunks...")
                chunk_files = split_fastq_files(r1_file, r2_file, num_chunks, chunks_dir)
                print(f"✓ Created {len(chunk_files)} chunks")
                
                # Process chunks sequentially
                chunk_bams = []
                for idx, (chunk_r1, chunk_r2) in enumerate(chunk_files):
                    chunk_bam = os.path.join(chunks_dir, f"chunk_{idx}.bam")
                    print(f"\n🔄 Processing chunk {idx + 1}/{len(chunk_files)}...")
                    
                    process_chunk(
                        chunk_r1,
                        chunk_r2,
                        chunk_bam,
                        idx0_file,
                        idx_mk_file,
                        ref_file,
                        fwd_lib,
                        max_mismatches,
                        threads,
                        min_alignment_length,
                        min_mapping_ratio,
                        idx,
                    )
                    chunk_bams.append(chunk_bam)
                    print(f"✓ Chunk {idx + 1} completed")
                
                # Merge all chunk BAMs into final output
                print(f"\n🔗 Merging {len(chunk_bams)} BAM files...")
                pysam.cat("-o", output_file, *chunk_bams)
                print(f"✓ Merged BAM saved to: {output_file}")
                
            finally:
                # Clean up chunks directory
                if os.path.exists(chunks_dir):
                    shutil.rmtree(chunks_dir)
            
            return

        # Create BAM header
        header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": []}
        for name, seq, *_ in mp.fastx_read(ref_file):
            header["SQ"].append({"SN": name, "LN": len(seq)})

        # Write BAM file (single chunk processing)
        with pysam.AlignmentFile(output_file, "wb", header=header) as bam_out:
            if r2_file is None:
                # Single-end
                with Progress() as progress:
                    task = progress.add_task("Mapping reads: 0 (0:00:00)", total=None)
                    count = 0
                    for name1, seq1, qua1 in mp.fastx_read(r1_file):
                        # Get mapped reads
                        mapped = run_mapping(
                            name1,
                            seq1,
                            None,
                            qua1,
                            None,
                            idx0,
                            idx_mk,
                            fwd_lib,
                            max_mismatches,
                            min_alignment_length,
                            min_mapping_ratio,
                        )

                        # Write to BAM
                        for i, item in enumerate(mapped):
                            map1 = item[1]
                            a1 = create_bam_record(
                                bam_out.header, map1, is_secondary=(i > 0)
                            )
                            bam_out.write(a1)

                        count += 1
                        if count % 10 == 0:
                            progress.update(
                                task,
                                description=f"Mapping reads: {count} ({format_duration(progress.tasks[task].elapsed)})",
                            )
            else:
                # Paired-end
                with Progress() as progress:
                    task = progress.add_task("Mapping reads: 0 (0:00:00)", total=None)
                    count = 0
                    for (name1, seq1, qua1), (name2, seq2, qua2) in zip(
                        mp.fastx_read(r1_file), mp.fastx_read(r2_file)
                    ):
                        if name1 != name2:
                            raise ValueError("r1 and r2 not in the same order")

                        # Get mapped reads
                        mapped = run_mapping(
                            name1,
                            seq1,
                            seq2,
                            qua1,
                            qua2,
                            idx0,
                            idx_mk,
                            fwd_lib,
                            max_mismatches,
                            min_alignment_length,
                            min_mapping_ratio,
                        )

                        # Write to BAM
                        for i, item in enumerate(mapped):
                            map1, map2 = item[1], item[2]
                            a1 = create_bam_record(
                                bam_out.header, map1, is_secondary=(i > 0)
                            )
                            a2 = create_bam_record(
                                bam_out.header, map2, is_secondary=(i > 0)
                            )
                            bam_out.write(a1)
                            bam_out.write(a2)

                        count += 1
                        if count % 10 == 0:
                            progress.update(
                                task,
                                description=f"Mapping reads: {count:,} ({format_duration(progress.tasks[task].elapsed)})",
                            )
    finally:
        # Clean up temporary directory only if we created it
        if use_temp_dir and os.path.exists(index_base_dir):
            shutil.rmtree(index_base_dir)
