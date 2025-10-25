#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-01-30 15:55

"""convert A->G, C->T in DNA sequence."""

import dnaio
from . import seqops


def convert_file(
    input_file: str,
    output_file: str,
    base_from: str,
    base_to: str,
    include_ys_tag: bool = True,
):
    """
    Convert DNA sequence from base_from to base_to.

    Args:
        input_file: Path to input FASTA/FASTQ file
        output_file: Path to output converted file
        base_from: Bases to convert from (e.g., "AC" for MK conversion)
        base_to: Bases to convert to (e.g., "GT" for MK conversion)
        include_ys_tag: Whether to include YS:Z tag with original sequence
    
    Examples:
        MK conversion (A->G, C->T): convert_file(in, out, "AC", "GT")
        KM conversion (G->A, T->C): convert_file(in, out, "GT", "AC")
    """
    with dnaio.open(input_file, mode="r") as fi, dnaio.open(output_file, mode="w") as fo:
        for read in fi:
            n = read.name.split()[0]

            if include_ys_tag:
                n = f"{n} YS:Z:{read.sequence}"

            r = dnaio.SequenceRecord(
                name=n,
                sequence=seqops.fast_base_conversion(read.sequence, base_from, base_to),
                qualities=read.qualities,
            )
            fo.write(r)


def convert_file_dual(
    input_file: str,
    output_X2Y_file: str,
    output_Y2X_file: str,
    base_from: str,
    base_to: str,
    include_ys_tag: bool = True,
):
    """
    Convert DNA sequence from base_from to base_to and vice versa (dual conversion).

    Args:
        input_file: Path to input FASTA/FASTQ file
        output_X2Y_file: Path to output file for base_from -> base_to conversion
        output_Y2X_file: Path to output file for base_to -> base_from conversion
        base_from: Bases to convert from (e.g., "AC")
        base_to: Bases to convert to (e.g., "GT")
        include_ys_tag: Whether to include YS:Z tag with original sequence
    
    Examples:
        base_from: "ACGT"
        base_to: "GTGT"
    """
    with dnaio.open(input_file, mode="r") as fi, dnaio.open(
        output_X2Y_file, mode="w"
    ) as fo_x2y, dnaio.open(output_Y2X_file, mode="w") as fo_y2x:
        for read in fi:
            n = read.name.split()[0]

            if include_ys_tag:
                n_x2y = f"{n} YS:Z:{read.sequence}"
                n_y2x = f"{n} YS:Z:{read.sequence}"
            else:
                n_x2y = n
                n_y2x = n

            r_x2y = dnaio.SequenceRecord(
                name=n_x2y,
                sequence=seqops.fast_base_conversion(read.sequence, base_from, base_to),
                qualities=read.qualities,
            )
            fo_x2y.write(r_x2y)

            r_y2x = dnaio.SequenceRecord(
                name=n_y2x,
                sequence=seqops.fast_base_conversion(read.sequence, base_to, base_from),
                qualities=read.qualities,
            )
            fo_y2x.write(r_y2x)


def mk_convert_file(
    input_file: str,
    output_M2K_file: str,
    output_K2M_file: str,
    include_ys_tag: bool = True,
):
    """Legacy function for dual MK/KM conversion."""
    convert_file_dual(
        input_file,
        output_M2K_file,
        output_K2M_file,
        "AC",
        "GT",
        include_ys_tag=include_ys_tag,
    )
