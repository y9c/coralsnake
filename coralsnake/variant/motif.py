#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Migrated from the standalone `variant` package (`variant motif`).
# pyfaidx is replaced with pysam.FastaFile (already a coralsnake dependency);
# gzip I/O and the fast reverse-complement are reused from coralsnake.

import sys

import pysam

from ..utils import get_logger, reverse_complement

LOGGER = get_logger(__name__)


def get_motif(fasta: pysam.FastaFile, chrom, chrom_len, pos, strand, lpad, rpad):
    """Fetch ``lpad + 1 + rpad`` bases centred on ``pos`` (1-based coordinate).

    Strand-aware: minus-strand sites are reverse-complemented. Out-of-bound
    positions are padded with ``N``.
    """
    if not str(pos).isdigit():
        LOGGER.error(f"Position {pos} is not a number!")
        sys.exit(1)
    if strand not in ["+", "-"]:
        LOGGER.error(f"Strand {strand} is not + or -!")
        sys.exit(1)

    # pos is 1-based, convert to 0-based
    pos = int(pos) - 1

    # Determine the (0-based, half-open) window, padding overhangs with N.
    if pos - lpad >= 0 and pos + rpad < chrom_len:
        start = pos - lpad
        end = pos + rpad + 1
        lfill = 0
        rfill = 0
    elif pos - lpad < 0 and pos + rpad < chrom_len:
        start = 0
        end = pos + rpad + 1
        lfill = lpad - pos
        rfill = 0
    elif pos - lpad >= 0 and pos + rpad >= chrom_len:
        start = pos - lpad
        end = chrom_len
        lfill = 0
        rfill = rpad - (chrom_len - pos)
    else:
        start = 0
        end = chrom_len
        lfill = lpad
        rfill = rpad - (chrom_len - pos)

    seq = fasta.fetch(chrom, start, end)
    if strand == "+":
        sequence = "N" * lfill + seq + "N" * rfill
    else:
        sequence = "N" * rfill + reverse_complement(seq) + "N" * lfill

    return sequence


def run_motif(
    input_file,
    output_file,
    fasta_path,
    lpad,
    rpad,
    with_header,
    columns,
    to_upper=True,
    wrap_site=True,
):
    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    columns_index_mapper = dict(zip(["chrom", "pos", "strand"], columns_index))
    strandness = "strand" in columns_index_mapper

    import gzip

    def open_in(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

    def open_out(path):
        return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

    with (
        open_in(input_file) as input_handle,
        open_out(output_file) as output_handle,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        chrom_len_mapper = dict(zip(fasta.references, fasta.lengths))

        def parse_line(input_cols):
            chrom_name = input_cols[columns_index_mapper["chrom"]]
            chrom_len = chrom_len_mapper[chrom_name]
            strand = input_cols[columns_index_mapper["strand"]] if strandness else "+"
            m = get_motif(
                fasta,
                chrom_name,
                chrom_len,
                input_cols[columns_index_mapper["pos"]],
                strand,
                lpad,
                rpad,
            )
            if to_upper:
                m = m.upper()
            if wrap_site:
                if strand == "+":
                    m = m[:lpad] + "[" + m[lpad] + "]" + m[lpad + 1 :]
                else:
                    m = m[:rpad] + "[" + m[rpad] + "]" + m[rpad + 1 :]

            output_handle.write(col_sep.join(input_cols + [m]) + "\n")

        # read first line and check column count
        input_cols = input_handle.readline().strip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(input_cols) - 1:
            LOGGER.error(f"Input file only have {len(input_cols)} columns!")
            sys.exit(1)
        if with_header:
            input_header = input_cols
        else:
            input_header = ["."] * len(input_cols)
            for n, i in columns_index_mapper.items():
                input_header[i] = n
        header_line = col_sep.join(input_header + ["motif"]) + "\n"
        # output header column only if input file is with header
        if with_header:
            output_handle.write(header_line)

        if not with_header:
            parse_line(input_cols)
        for line in input_handle:
            input_cols = line.strip("\n").split(col_sep)
            parse_line(input_cols)
