#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-01-30 15:55

"""convert A->G, C->T in DNA sequence."""

import sys

import dnaio

MK_BASE_MAPPER = str.maketrans("ACGTacgt", "GTGTgtgt")
KM_BASE_MAPPER = str.maketrans("ACGTacgt", "ACACacac")


def mk_conversion(seq):
    return seq.translate(MK_BASE_MAPPER)


def km_conversion(seq):
    return seq.translate(KM_BASE_MAPPER)


def convert_file(input_file, output_MK_file, output_KM_file):
    with dnaio.open(input_file, mode="r") as fi, dnaio.open(
        output_MK_file, mode="w"
    ) as fo_mk, dnaio.open(output_KM_file, mode="w") as fo_km:
        for read in fi:
            n = read.name.split()[0]
            n_mk = f"{n} YS:Z:{read.sequence}\tST:A:+"
            r_mk = dnaio.SequenceRecord(
                name=n_mk,
                sequence=mk_conversion(read.sequence),
                qualities=read.qualities,
            )
            fo_mk.write(r_mk)

            n_km = f"{n} YS:Z:{read.sequence}\tST:A:-"
            r_km = dnaio.SequenceRecord(
                name=n_km,
                sequence=km_conversion(read.sequence),
                qualities=read.qualities,
            )
            fo_km.write(r_km)
