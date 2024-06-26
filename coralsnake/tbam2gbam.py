#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-23 18:01


import bisect
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
    if flag & 16:
        return (flag & ~16) | 32
    elif flag & 32:
        return (flag & ~32) | 16
    else:
        return flag ^ 16


@lru_cache(maxsize=10000)
def reverse_md(md: str) -> str:
    parts = []
    num = ""
    for char in md:
        if char.isdigit():
            num += char
        else:
            if num:
                parts.append(num)
                num = ""
            if char == "^":
                deletion = "^"
                while md and md[0].isalpha():
                    deletion += md[0]
                    md = md[1:]
                parts.append(deletion[::-1])
            else:
                parts.append(char.translate(COMP))
    if num:
        parts.append(num)
    return "".join(reversed(parts))


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
    align: pysam.AlignedSegment, transcript: Transcript, header: pysam.AlignmentHeader
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
    annot: dict[str, dict[str, Transcript]],
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

    transcript = next(iter(annot[align.reference_name].values()))
    new_align = remap_to_genome(align, transcript, genome_header)
    return new_align


def convert_bam(
    input_bam: str,
    output_bam: str,
    annotation_file: str,
    faidx_file: str,
    threads: int = 8,
    sort: bool = False,
):
    LOGGER.info("Loading annotation...")
    annot = load_annotation(annotation_file)
    LOGGER.info("Loading fasta index...")
    faidx = load_faidx(faidx_file)

    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        new_header = in_bam.header.to_dict()
        # mark as unsorted
        new_header["HD"]["SO"] = "unsorted"
        new_header["SQ"] = [
            {"SN": chrom, "LN": length} for chrom, length in faidx.items()
        ]
        genome_header = pysam.AlignmentHeader.from_dict(new_header)

        with pysam.AlignmentFile(
            output_bam,
            "wb" if output_bam.endswith(".bam") else "w",
            header=genome_header,
        ) as out_bam:
            for align in track(in_bam, description="Processing..."):
                new_align = parse_alignment(align, annot, genome_header)
                out_bam.write(new_align)

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
