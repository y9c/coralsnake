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
from ncls import NCLS32
from rich.progress import track
from xopen import xopen


# Function to parse exon data into a structured dictionary
def parse_annot(tx_file, cache=True):
    # check if pickle file exists, load it if it does
    if cache and os.path.exists(tx_file + ".pickle"):
        with open(tx_file + ".pickle", "rb") as f:
            data = pickle.load(f)
            exon_tree_by_chrom_strand = data["tree"]
            for (chromosome, strand), (
                starts,
                ends,
                rids,
            ) in exon_tree_by_chrom_strand.items():
                exon_tree_by_chrom_strand[(chromosome, strand)] = NCLS32(
                    starts, ends, rids
                )
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
            info[rid] = (gene_id, transcript_id)
            if (chromosome, strand) not in exons_by_chrom_strand:
                exons_by_chrom_strand[(chromosome, strand)] = []
            for exon_range in exon_positions.split(","):
                start, end = map(int, exon_range.split("-"))
                exons_by_chrom_strand[(chromosome, strand)].append((start, end, rid))
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
            pickle.dump(
                {"tree": exon_tree_by_chrom_strand, "info": info},
                open(tx_file + ".pickle", "wb"),
            )

    for (chromosome, strand), (starts, ends, rids) in exon_tree_by_chrom_strand.items():
        exon_tree_by_chrom_strand[(chromosome, strand)] = NCLS32(starts, ends, rids)
    return exon_tree_by_chrom_strand, info


def annot_file(
    input_file, output_file, annot_file, cols=None, keep_na=True, skip_header=False
):
    tree_by_chrom_strand, info = parse_annot(annot_file)
    if cols is None:
        cols = [0, 1, 2]
    else:
        cols = [int(i) - 1 for i in cols.split(",")]
    with xopen(input_file, "rt") as fi, xopen(output_file, "wt") as fo:
        if skip_header:
            next(fi)
        for line in track(fi, description="Processing sites"):
            line = line.strip("\n")
            records = line.split("\t")
            chromosome, position, strand = [records[i] for i in cols]
            position = int(position) - 1
            tree = tree_by_chrom_strand.get((chromosome, strand))
            founded = False
            if tree:
                for _, _, rid in tree.find_overlap(position, position + 1):
                    if info.get(rid):
                        founded = True
                        gene_id, transcript_id = info[rid]
                        fo.write(f"{line}\t{gene_id}\t{transcript_id}\n")
            if not founded and keep_na:
                fo.write(f"{line}\t.\t.\n")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("input_file", help="site file")
    ap.add_argument("output_file", help="output file")
    ap.add_argument("annot_file", help="annotation file")
    args = ap.parse_args()
    annot_file(args.input_file, args.output_file, args.annot_file)
