#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-01-30 15:55

"""convert A->G, C->T in DNA sequence."""

import dnaio


def create_base_mapper(base_from: str, base_to: str) -> dict[int, int]:
    if len(base_from) != len(base_to):
        raise ValueError("base_from and base_to must be of the same length")

    # Create mapping for both lower and upper case
    base_from_lower = base_from.lower()
    base_from_upper = base_from.upper()
    base_to_lower = base_to.lower()
    base_to_upper = base_to.upper()

    full_base_from = base_from_lower + base_from_upper
    full_base_to = base_to_lower + base_to_upper

    return str.maketrans(full_base_from, full_base_to)


def base_conversion(seq: str, base_mapper: dict[int, int]) -> str:
    return seq.translate(base_mapper)


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
    X2Y_BASE_MAPPER = create_base_mapper(base_from, base_to)
    Y2X_BASE_MAPPER = create_base_mapper(base_to, base_from)

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
                sequence=base_conversion(read.sequence, X2Y_BASE_MAPPER),
                qualities=read.qualities,
            )
            fo_x2y.write(r_x2y)

            r_y2x = dnaio.SequenceRecord(
                name=n_y2x,
                sequence=base_conversion(read.sequence, Y2X_BASE_MAPPER),
                qualities=read.qualities,
            )
            fo_y2x.write(r_y2x)


# M to K (AC to GT)
# M (aMino)	A/C
# K (Keto)	G/T


def mk_conversion(seq: str) -> str:
    return base_conversion(seq, create_base_mapper("AC", "GT"))


def km_conversion(seq: str) -> str:
    return base_conversion(seq, create_base_mapper("GT", "AC"))


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
