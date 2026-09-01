#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-29 19:54


import os
import pickle

import numpy as np
from ruranges.numpy import overlaps
from xopen import xopen


# Function to parse exon data into a structured dictionary
def parse_annot_file(tx_file, cache):
    # check if pickle file exists, load it if it does (but never reuse a cache
    # that is older than the table itself)
    pickle_file = tx_file + ".pickle"
    if (
        cache
        and os.path.exists(pickle_file)
        and os.path.getmtime(pickle_file) >= os.path.getmtime(tx_file)
    ):
        with open(pickle_file, "rb") as f:
            data = pickle.load(f)
            exon_tree_by_chrom_strand = data["tree"]
            return exon_tree_by_chrom_strand, data["info"]

    exons_by_chrom_strand = {}
    info = {}
    with open(tx_file, "r") as f:
        names = f.readline().strip("\n").split("\t")
        chrom_idx = names.index("chrom")
        strand_idx = names.index("strand")
        spans_idx = names.index("spans")
        gene_idx = names.index("gene_id")
        transcript_idx = names.index("transcript_id")

        rid = 0
        for line in f:
            record = line.strip("\n").split("\t")
            chromosome = record[chrom_idx]
            strand = record[strand_idx]
            exon_positions = record[spans_idx]
            gene_id = record[gene_idx]
            transcript_id = record[transcript_idx]
            exon_shift = 0
            if (chromosome, strand) not in exons_by_chrom_strand:
                exons_by_chrom_strand[(chromosome, strand)] = []
            for exon_range in exon_positions.split(","):
                start, end = map(int, exon_range.split("-"))
                # annotation file is 1-based, convert to 0-based
                start -= 1
                exons_by_chrom_strand[(chromosome, strand)].append((start, end, rid))
                info[rid] = (gene_id, transcript_id, exon_shift)
                exon_shift += end - start
                rid += 1

    exon_tree_by_chrom_strand = {}
    for (chromosome, strand), exons in exons_by_chrom_strand.items():
        starts, ends, rid = list(zip(*exons))
        starts = np.int32(starts)
        ends = np.int32(ends)
        rids = np.int64(rid)
        exon_tree_by_chrom_strand[(chromosome, strand)] = (starts, ends, rids)
    # save the parsed data to a pickle file
    if cache:
        with open(tx_file + ".pickle", "wb") as f:
            pickle.dump({"tree": exon_tree_by_chrom_strand, "info": info}, f)

    return exon_tree_by_chrom_strand, info


def _read_sites(fi, cols, skip_header):
    """Read the input site file into aligned column lists (streaming-friendly)."""
    chroms = []
    positions = []
    strands = []
    lines = []
    if skip_header:
        next(fi)
    for raw in fi:
        line = raw.strip("\n")
        records = line.split("\t")
        try:
            chrom = records[cols[0]]
            pos = int(records[cols[1]]) - 1
            strand = records[cols[2]]
        except (ValueError, IndexError):
            continue  # skip a malformed row instead of aborting the run
        chroms.append(chrom)
        positions.append(pos)
        strands.append(strand)
        lines.append(line)
    return lines, chroms, positions, strands


def _annotate_batch(
    lines, chroms, positions, strands, tree_by_chrom_strand, info, batch_size=200000
):
    """Annotate sites in vectorized batches instead of one overlaps() call each.

    Returns a list (parallel to ``lines``) of records: each is a list of
    (gene_id, transcript_id, transcript_pos) tuples in the order returned by
    the overlap engine (matching per-site behavior).
    """
    n = len(lines)
    results = [[] for _ in range(n)]

    # Group (chrom, strand) -> stable group id and the concatenated exon arrays.
    # Only sites whose (chrom, strand) has an exon tree can be annotated.
    tree_keys = list(tree_by_chrom_strand.keys())
    tree_id = {k: i for i, k in enumerate(tree_keys)}

    exon_starts = []
    exon_ends = []
    exon_rids = []
    exon_groups = []
    for key in tree_keys:
        starts, ends, rids = tree_by_chrom_strand[key]
        g = tree_id[key]
        exon_starts.append(starts)
        exon_ends.append(ends)
        exon_rids.append(rids)
        exon_groups.append(np.full(len(starts), g, dtype=np.uint32))
    if not exon_starts:
        # No annotations at all -> every site is unannotated.
        return results

    ex_starts = np.concatenate(exon_starts).astype(np.int32)
    ex_ends = np.concatenate(exon_ends).astype(np.int32)
    ex_rids = np.concatenate(exon_rids).astype(np.int64)
    ex_groups = np.concatenate(exon_groups)

    for start in range(0, n, batch_size):
        stop = min(start + batch_size, n)
        # Sites in this batch that have an exon tree.
        q_idx = [i for i in range(start, stop) if (chroms[i], strands[i]) in tree_id]
        if not q_idx:
            continue
        q_starts = np.array([positions[i] for i in q_idx], dtype=np.int32)
        q_ends = q_starts + 1
        q_groups = np.array(
            [tree_id[(chroms[i], strands[i])] for i in q_idx], dtype=np.uint32
        )

        idx_site, idx_exon = overlaps(
            starts=q_starts,
            ends=q_ends,
            starts2=ex_starts,
            ends2=ex_ends,
            groups=q_groups,
            groups2=ex_groups,
        )
        # Group matches by their batch position (idx_site indexes into q_idx).
        for k in range(len(idx_site)):
            site_i = q_idx[int(idx_site[k])]
            j = int(idx_exon[k])
            rid = int(ex_rids[j])
            meta = info.get(rid)
            if not meta:
                continue
            gene_id, transcript_id, exon_shift = meta
            exon_start = int(ex_starts[j])
            exon_end = int(ex_ends[j])
            position = positions[site_i]
            if strands[site_i] == "+":
                transcript_pos = position - exon_start + exon_shift
            else:
                transcript_pos = exon_end - 1 - position + exon_shift
            results[site_i].append((gene_id, transcript_id, transcript_pos))
    return results


def run_annot(
    input_file,
    output_file,
    annot_file,
    cols=None,
    keep_na=True,
    collapse_annot=False,
    add_count=False,
    skip_header=False,
):
    cache = True
    # if not collapse_annot, we can add a column of annotation count
    if collapse_annot and add_count:
        raise ValueError("--collapse-annot and --add-count cannot be both True")
    tree_by_chrom_strand, info = parse_annot_file(annot_file, cache)
    if cols is None:
        cols = [0, 1, 2]
    else:
        cols = [int(i) - 1 for i in cols.split(",")]
    with xopen(input_file, "rt") as fi, xopen(output_file, "wt") as fo:
        lines, chroms, positions, strands = _read_sites(fi, cols, skip_header)
        results = _annotate_batch(
            lines, chroms, positions, strands, tree_by_chrom_strand, info
        )
        for line, annot_list in zip(lines, results):
            if annot_list:
                if collapse_annot:
                    gene_id_join = ",".join([x[0] for x in annot_list])
                    transcript_id_join = ",".join([x[1] for x in annot_list])
                    transcript_pos_join = ",".join([str(x[2]) for x in annot_list])
                    fo.write(
                        f"{line}\t{gene_id_join}\t{transcript_id_join}\t{transcript_pos_join}\n"
                    )
                else:
                    for gene_id, transcript_id, transcript_pos in annot_list:
                        fo.write(
                            f"{line}\t{gene_id}\t{transcript_id}\t{transcript_pos}\n"
                            if not add_count
                            else f"{line}\t{gene_id}\t{transcript_id}\t{transcript_pos}\t{len(annot_list)}\n"
                        )
            else:
                if keep_na:
                    fo.write(
                        f"{line}\t.\t.\t.\n"
                        if not add_count
                        else f"{line}\t.\t.\t.\t0\n"
                    )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("input_file", help="site file")
    ap.add_argument("output_file", help="output file")
    ap.add_argument("annot_file", help="annotation file")
    args = ap.parse_args()
    run_annot(args.input_file, args.output_file, args.annot_file)
