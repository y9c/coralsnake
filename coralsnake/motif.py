#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Genomic motif fetch, fused into coralsnake from the standalone `variant`
# package (`variant motif`). pyfaidx is replaced by pysam.FastaFile (already a
# coralsnake dependency); gzip/stdin-stdout I/O uses xopen.

import pysam

from .utils import get_logger, reverse_complement

LOGGER = get_logger(__name__)


def get_motif(fasta: pysam.FastaFile, chrom, chrom_len, pos, strand, lpad, rpad):
    """Fetch ``lpad + 1 + rpad`` bases centred on ``pos`` (1-based coordinate).

    Strand-aware: minus-strand sites are reverse-complemented. Out-of-bound
    positions are padded with ``N``.
    """
    if not str(pos).isdigit():
        raise ValueError(f"Position {pos} is not a number!")
    if strand not in ["+", "-"]:
        raise ValueError(f"Strand {strand} is not + or -!")

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


def _wrap_site(motif, strand, lpad, rpad):
    """Wrap the motif centre site in ``[...]``."""
    if strand == "+":
        return motif[:lpad] + "[" + motif[lpad] + "]" + motif[lpad + 1 :]
    return motif[:rpad] + "[" + motif[rpad] + "]" + motif[rpad + 1 :]


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
    """Fetch a motif for every variant site; write to ``output_file``."""
    from xopen import xopen

    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    columns_index_mapper = dict(zip(["chrom", "pos", "strand"], columns_index))
    strand_col = columns_index_mapper.get("strand")

    with (
        xopen(input_file, "rt") as input_handle,
        xopen(output_file, "wt") as output_handle,
        pysam.FastaFile(fasta_path) as fasta,
    ):
        chrom_len_mapper = dict(zip(fasta.references, fasta.lengths))

        def parse_line(input_cols):
            chrom_name = input_cols[columns_index_mapper["chrom"]]
            chrom_len = chrom_len_mapper[chrom_name]
            strand = input_cols[strand_col] if strand_col is not None else "+"
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
                m = _wrap_site(m, strand, lpad, rpad)
            output_handle.write(col_sep.join(input_cols + [m]) + "\n")

        # Read the first line to infer header / column layout.
        first = input_handle.readline().rstrip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(first) - 1:
            raise ValueError(f"Input file only has {len(first)} columns!")
        if with_header:
            input_header = first
            output_handle.write(col_sep.join(input_header + ["motif"]) + "\n")
        else:
            parse_line(first)

        for line in input_handle:
            if not line.strip():
                continue
            parse_line(line.rstrip("\n").split(col_sep))
