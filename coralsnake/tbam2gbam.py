#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-23 18:01


import bisect
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache

import pysam
from rich.progress import track

from .utils import (
    Transcript,
    get_logger,
    load_annotation,
    load_faidx,
    reverse_complement,
)

LOGGER = get_logger(__name__)

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@lru_cache(maxsize=10000)
def flip_flag(flag):
    if flag & 1:  # paired: toggle both read-reverse (16) and mate-reverse (32)
        return flag ^ 0x30
    else:  # unpaired: toggle only read-reverse (16)
        return flag ^ 16


@lru_cache(maxsize=10000)
def reverse_md(md: str) -> str:
    # Parse MD into components: numbers, mismatches, deletions
    parts = []
    import re

    # Match numbers, mismatches (single base), or deletions (^ followed by bases)
    for m in re.finditer(r"([0-9]+)|([A-Z])|(\^[A-Z]+)", md):
        if m.group(1) is not None:  # Number
            parts.append(int(m.group(1)))
        elif m.group(2):  # Mismatch base
            parts.append(m.group(2).translate(COMP))
        elif m.group(3):  # Deletion
            # Complement the deletion bases and reverse them
            del_bases = m.group(3)[1:]
            parts.append("^" + del_bases.translate(COMP)[::-1])

    # Reverse the parts list
    parts.reverse()

    # Merge adjacent numbers, ensuring 0s are kept between non-numbers
    merged_parts = []
    current_num = 0
    for p in parts:
        if isinstance(p, int):
            current_num += p
        else:
            merged_parts.append(str(current_num))
            current_num = 0
            merged_parts.append(p)
    merged_parts.append(str(current_num))
    return "".join(merged_parts)


@lru_cache(maxsize=100000)
def transcript_to_genome(
    transcript_pos: int, transcript: Transcript
) -> tuple[int, int]:
    """
    transcript_pos is 0-based
    """
    if transcript_pos > transcript.length or transcript_pos < 0:
        raise ValueError("Transcript position is out of range")
    # bisect_right dose not include the right bound
    # in math notation, it is [a, b)
    exon_index = bisect.bisect_right(transcript.cum_exon_lens, transcript_pos)
    if exon_index == 0:
        return transcript.exons_forwards[0].start + transcript_pos, 0
    offset = transcript_pos - transcript.cum_exon_lens[exon_index - 1]
    return transcript.exons_forwards[exon_index].start + offset, exon_index


def remap_to_genome(
    align: pysam.AlignedSegment,
    header: pysam.AlignmentHeader,
    transcript: Transcript,
    next_transcript: Transcript,
) -> pysam.AlignedSegment:
    new_align = pysam.AlignedSegment(header=header)
    new_align.query_name = align.query_name
    new_align.flag = align.flag

    if align.cigartuples is None:
        raise ValueError("CIGAR string is required")

    for tag in align.get_tags():
        new_align.set_tag(*tag)
    if transcript.strand == "-":
        new_align.flag = flip_flag(align.flag)
        new_align.query_sequence = reverse_complement(align.query_sequence)
        new_align.query_qualities = (
            align.query_qualities[::-1] if align.query_qualities else None
        )
        cigartuples = align.cigartuples[::-1]
        if new_align.has_tag("MD"):
            new_align.set_tag("MD", reverse_md(new_align.get_tag("MD")))
        transcript_pos = transcript.length - align.reference_end
    else:
        new_align.query_sequence = align.query_sequence
        new_align.query_qualities = align.query_qualities
        cigartuples = align.cigartuples
        transcript_pos = align.reference_start

    genome_pos, exon_index = transcript_to_genome(transcript_pos, transcript)
    new_align.reference_name = transcript.chrom
    if next_transcript is not None:
        new_align.next_reference_name = next_transcript.chrom
    new_align.reference_start = genome_pos

    new_cigar = []
    for op, length in cigartuples:
        if op in (0, 2):  # Match, deletion
            while length > 0:
                current_exon = transcript.exons_forwards[exon_index]
                exon_remaining = current_exon.end - genome_pos
                if length <= exon_remaining:
                    new_cigar.append((op, length))
                    genome_pos += length
                    length = 0
                else:
                    if exon_remaining > 0:
                        new_cigar.append((op, exon_remaining))
                        genome_pos += exon_remaining
                        length -= exon_remaining
                    if exon_index + 1 < len(transcript.exons_forwards):
                        next_exon = transcript.exons_forwards[exon_index + 1]
                        intron_length = next_exon.start - current_exon.end
                        new_cigar.append((3, intron_length))
                        genome_pos = next_exon.start
                        exon_index += 1
                    else:
                        break
        else:  # Insertion or soft-clip
            new_cigar.append((op, length))

    new_align.cigartuples = new_cigar
    new_align.mapping_quality = align.mapping_quality

    return new_align


def parse_alignment(
    align: pysam.AlignedSegment,
    annot: dict[str, Transcript],
    genome_header: pysam.AlignmentHeader,
) -> pysam.AlignedSegment:
    if (
        align.is_unmapped
        or align.reference_name is None
        or align.reference_name not in annot
    ):
        align.reference_id = -1
        align.reference_start = -1
        align.next_reference_id = -1
        align.next_reference_start = -1
        return align

    transcript = annot[align.reference_name]
    if align.next_reference_name in annot:
        next_transcript = annot[align.next_reference_name]
    else:
        next_transcript = None
    new_align = remap_to_genome(align, genome_header, transcript, next_transcript)
    return new_align


# Worker global state for multiprocessing
_WORKER_ANNOT = None
_WORKER_GENOME_HEADER = None
_WORKER_TRANSCRIPT_HEADER = None


def _init_worker(annotation_file, genome_header_dict, transcript_header_dict):
    """Initialize worker with shared annotation and headers."""
    global _WORKER_ANNOT, _WORKER_GENOME_HEADER, _WORKER_TRANSCRIPT_HEADER
    raw_annot = load_annotation(annotation_file)
    # Build flat mapping for transcript lookups
    _WORKER_ANNOT = {}
    for g_id in raw_annot:
        for t_id, tx in raw_annot[g_id].items():
            _WORKER_ANNOT[t_id] = tx

    _WORKER_GENOME_HEADER = pysam.AlignmentHeader.from_dict(genome_header_dict)
    _WORKER_TRANSCRIPT_HEADER = pysam.AlignmentHeader.from_dict(transcript_header_dict)


def _process_chunk(align_strings):
    """Process a chunk of alignments (serialized as strings)."""
    results = []
    for s in align_strings:
        # Deserialize from string using transcript header
        align = pysam.AlignedSegment.fromstring(s, _WORKER_TRANSCRIPT_HEADER)
        # Remap to genome
        new_align = parse_alignment(align, _WORKER_ANNOT, _WORKER_GENOME_HEADER)
        # Serialize result to string for return
        results.append(new_align.to_string())
    return results


def convert_bam(
    input_bam: str,
    output_bam: str,
    annotation_file: str,
    faidx_file: str,
    threads: int = 8,
    sort: bool = False,
):
    LOGGER.info("Loading fasta index...")
    faidx = load_faidx(faidx_file)

    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        transcript_header_dict = in_bam.header.to_dict()
        new_header = in_bam.header.to_dict()
        if "HD" not in new_header:
            new_header["HD"] = {"VN": "1.4"}
        # mark as unsorted
        new_header["HD"]["SO"] = "unsorted"
        new_header["SQ"] = [
            {"SN": chrom, "LN": length} for chrom, length in faidx.items()
        ]
        genome_header_dict = new_header
        genome_header = pysam.AlignmentHeader.from_dict(genome_header_dict)

        with pysam.AlignmentFile(
            output_bam,
            "wb" if output_bam.endswith(".bam") else "w",
            header=genome_header,
        ) as out_bam:
            if threads <= 1:
                LOGGER.info("Loading annotation for single-threaded processing...")
                raw_annot = load_annotation(annotation_file)
                annot = {}
                for g_id in raw_annot:
                    for t_id, tx in raw_annot[g_id].items():
                        annot[t_id] = tx
                for align in track(in_bam, description="Processing..."):
                    new_align = parse_alignment(align, annot, genome_header)
                    out_bam.write(new_align)
            else:
                LOGGER.info(f"Using {threads} workers for parallel processing...")
                chunk_size = 5000
                # Keep a limited number of futures in flight to save memory
                max_queue = threads * 2
                futures = []
                chunk = []

                with ProcessPoolExecutor(
                    max_workers=threads,
                    mp_context=mp.get_context("spawn"),
                    initializer=_init_worker,
                    initargs=(
                        annotation_file,
                        genome_header_dict,
                        transcript_header_dict,
                    ),
                ) as executor:
                    for align in track(in_bam, description="Processing..."):
                        chunk.append(align.to_string())
                        if len(chunk) >= chunk_size:
                            futures.append(executor.submit(_process_chunk, list(chunk)))
                            chunk = []

                            # Drain results if queue is full
                            if len(futures) >= max_queue:
                                for res_str in futures.pop(0).result():
                                    res_align = pysam.AlignedSegment.fromstring(
                                        res_str, genome_header
                                    )
                                    out_bam.write(res_align)

                    # Submit final chunk
                    if chunk:
                        futures.append(executor.submit(_process_chunk, chunk))

                    # Drain remaining futures in order to maintain read sequence
                    for fut in futures:
                        for res_str in fut.result():
                            res_align = pysam.AlignedSegment.fromstring(
                                res_str, genome_header
                            )
                            out_bam.write(res_align)

    if sort:
        # skip sorting if output is not BAM file
        if not output_bam.endswith(".bam"):
            LOGGER.info("Output is not BAM file, skip sorting and indexing")
            return
        LOGGER.info("Sorting output BAM file...")
        import shutil
        import uuid

        internal_sorted_bam = f".{str(uuid.uuid4())[:8]}_{output_bam}"
        pysam.sort("-@", str(threads), "-o", internal_sorted_bam, output_bam)
        shutil.move(internal_sorted_bam, output_bam)
        pysam.index("-@", str(threads), output_bam)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Remap transcriptome BAM to genome BAM"
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file")
    parser.add_argument("-f", "--faidx", required=True, help="Fasta index file")
    # sort bam file
    parser.add_argument(
        "--sort", action="store_true", help="Sort output BAM file by coordinates"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=8, help="Number of threads"
    )
    args = parser.parse_args()

    convert_bam(
        args.input, args.output, args.annotation, args.faidx, args.threads, args.sort
    )
