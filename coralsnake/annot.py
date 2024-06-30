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
def parse_annot_file(tx_file, cache):
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
            pickle.dump(
                {"tree": exon_tree_by_chrom_strand, "info": info},
                open(tx_file + ".pickle", "wb"),
            )

    for (chromosome, strand), (starts, ends, rids) in exon_tree_by_chrom_strand.items():
        exon_tree_by_chrom_strand[(chromosome, strand)] = NCLS32(starts, ends, rids)
    return exon_tree_by_chrom_strand, info


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
    # TODO: add cache option
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
        if skip_header:
            next(fi)
        for line in track(fi, description="Processing sites"):
            line = line.strip("\n")
            records = line.split("\t")
            chromosome, position, strand = [records[i] for i in cols]
            position = int(position) - 1
            tree = tree_by_chrom_strand.get((chromosome, strand))
            annot_list = []
            if tree:
                for exon_start, exon_end, rid in tree.find_overlap(
                    position, position + 1
                ):
                    if info.get(rid):
                        gene_id, transcript_id, exon_shift = info[rid]
                        if strand == "+":
                            transcript_pos = position - exon_start + exon_shift
                        else:
                            transcript_pos = exon_end - position + exon_shift
                        annot_list.append((gene_id, transcript_id, transcript_pos))

            if len(annot_list) > 0:
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
