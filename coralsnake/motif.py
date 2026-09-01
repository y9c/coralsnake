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

# Sliding-window size for batched motif fetching (see ``_fetch_window``).
_WINDOW = 4_000_000


def _fetch_window(window_cache, fasta, chrom, chrom_len, start, end):
    """Return ``fasta[chrom][start:end]``, caching a larger surrounding window.

    ``run_motif`` processes each chromosome's sites in sorted order, so one
    fetch per ``_WINDOW`` bp replaces one faidx seek+fread per site (a large
    speedup when sites are spread across a big genome). ``window_cache`` is a
    ``{chrom: (w_start, w_end, seq)}`` dict holding at most one chromosome.
    """
    entry = window_cache.get(chrom)
    if entry is not None and entry[0] <= start and end <= entry[1]:
        w_start, _w_end, seq = entry
        return seq[start - w_start : end - w_start]
    w_end = min(chrom_len, max(end, start + _WINDOW))
    seq = fasta.fetch(chrom, start, w_end)
    if len(seq) < end - start:  # defensive: truncated fetch
        seq = seq + "N" * (end - start - len(seq))
    window_cache[chrom] = (start, w_end, seq)
    return seq[: end - start]


def get_motif(
    fasta: pysam.FastaFile,
    chrom,
    chrom_len,
    pos,
    strand,
    lpad,
    rpad,
    window_cache=None,
):
    """Fetch ``lpad + 1 + rpad`` bases centred on ``pos`` (1-based coordinate).

    Strand-aware: minus-strand sites are reverse-complemented. Out-of-bound
    positions are padded with ``N``. ``window_cache`` (optional) enables the
    batched-fetch path used by :func:`run_motif`.
    """
    if not str(pos).isdigit():
        raise ValueError(f"Position {pos} is not a number!")
    if strand not in ["+", "-"]:
        raise ValueError(f"Strand {strand} is not + or -!")

    # pos is 1-based, convert to 0-based
    pos = int(pos) - 1

    # Window centred on `pos` (0-based, half-open). Clamp to the chromosome and
    # pad any out-of-bound flank with N so the motif is always lpad+1+rpad long:
    #   left flank  = indices (pos - lpad) .. pos-1
    #   right flank = indices (pos + 1) .. (pos + rpad)
    start = max(0, pos - lpad)
    end = min(chrom_len, pos + rpad + 1)
    lfill = max(0, lpad - pos)            # leading bases that fall before 0
    rfill = max(0, rpad - (chrom_len - pos - 1))  # trailing bases past the end

    if start >= end:
        # every base of the window falls outside the contig -> all-N motif
        return "N" * (lpad + 1 + rpad)

    if window_cache is None:
        seq = fasta.fetch(chrom, start, end)
    else:
        seq = _fetch_window(window_cache, fasta, chrom, chrom_len, start, end)
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
    """Fetch a motif for every variant site; write to ``output_file``.

    Sites are collected, grouped by chromosome and processed in position order
    through a small sliding FASTA window (see ``_fetch_window``), so the many
    per-site faidx reads collapse into a few large sequential reads. Output is
    written back in the original input order. Validates every row up front so
    a malformed row aborts exactly where the sequential path would.
    """
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

        # Read the first line to infer header / column layout.
        first = input_handle.readline().rstrip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(first) - 1:
            raise ValueError(f"Input file only has {len(first)} columns!")
        if with_header:
            input_header = first
            output_handle.write(col_sep.join(input_header + ["motif"]) + "\n")

        # Collect + validate rows in input order (same errors as the
        # sequential path, raised at the same first-offending row).
        rows = []  # (input_cols, pos_int, strand)
        if not with_header:
            rows.append(first)
        for line in input_handle:
            if not line.strip():
                continue
            rows.append(line.rstrip("\n").split(col_sep))

        by_chrom = {}
        for r, cols in enumerate(rows):
            chrom_name = cols[columns_index_mapper["chrom"]]
            pos_str = cols[columns_index_mapper["pos"]]
            if not str(pos_str).isdigit():
                raise ValueError(f"Position {pos_str} is not a number!")
            strand = cols[strand_col] if strand_col is not None else "+"
            if strand not in ["+", "-"]:
                raise ValueError(f"Strand {strand} is not + or -!")
            by_chrom.setdefault(chrom_name, []).append(r)

        results = [None] * len(rows)
        window_cache = {}
        for chrom_name, rids in by_chrom.items():
            if chrom_name not in chrom_len_mapper:
                raise ValueError(
                    f"Chromosome {chrom_name!r} not found in the FASTA file"
                )
            chrom_len = chrom_len_mapper[chrom_name]
            # Position order maximises sliding-window reuse within one chrom.
            rids.sort(key=lambda r: int(rows[r][columns_index_mapper["pos"]]))
            for r in rids:
                cols = rows[r]
                strand = cols[strand_col] if strand_col is not None else "+"
                m = get_motif(
                    fasta,
                    chrom_name,
                    chrom_len,
                    cols[columns_index_mapper["pos"]],
                    strand,
                    lpad,
                    rpad,
                    window_cache=window_cache,
                )
                if to_upper:
                    m = m.upper()
                if wrap_site:
                    m = _wrap_site(m, strand, lpad, rpad)
                results[r] = m
            window_cache.pop(chrom_name, None)  # bound memory to one chrom

        for cols, m in zip(rows, results):
            output_handle.write(col_sep.join(cols + [m]) + "\n")
