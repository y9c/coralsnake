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
    output_X2Y_file: str,
    output_Y2X_file: str,
    base_from: str,
    base_to: str,
    include_ys_tag: bool = True,
):
    """
    Convert DNA sequence from base_from to base_to and vice versa.

    base_from: eg, "ACGT"
    base_to: eg, "GTGT"
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
    convert_file(
        input_file,
        output_M2K_file,
        output_K2M_file,
        "AC",
        "GT",
        include_ys_tag=include_ys_tag,
    )
