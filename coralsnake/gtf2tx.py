#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-05-01 18:23


import logging
from collections import defaultdict

from pyfaidx import Fasta
from rich.progress import track

logging.basicConfig(level=logging.INFO)


class TranscriptExon:
    def __init__(self, chrom=None, strand=None, gene_id=None, gene_name=None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chrom = chrom
        self.strand = strand
        self.span = {}
        self.priority = (10, 0)

    def add_span(self, gene_id, gene_name, chrom, strand, exon_number, start, end):
        self.gene_id = gene_id
        self.gene_name = gene_name
        if self.chrom is not None and self.chrom != chrom:
            raise ValueError("Chrom mismatch")
        self.chrom = chrom
        self.strand = strand
        self.span[exon_number] = (start, end)

    def calc_len(self):
        return sum([v[1] - v[0] + 1 for _, v in self.span.items()])

    def get_seq(self, fasta):
        seq = ""
        for _, v in sorted(self.span.items()):
            e = fasta[self.chrom][v[0] - 1 : v[1]]
            if self.strand == "-":
                e = e.reverse.complement
            seq += e.seq
        return seq.upper()

    def __repr__(self):
        # sort the exon number
        # span_str = ",".join([f"{v[0]}-{v[1]}" for _, v in sorted(self.span.items())])
        # do not sort the exon number
        span_str = ",".join([f"{v[0]}-{v[1]}" for v in self.span.values()])
        gene_name = "." if self.gene_name is None else self.gene_name
        return f"{gene_name}\t{self.chrom}\t{self.strand}\t{span_str}"


def parse_gtf_annot(annot):
    annot = annot.rstrip("\n").rstrip(";").split("; ")
    annot = [x.strip().split(" ", 1) for x in annot if x]
    # if the key is duplicate, join the values
    d = defaultdict(list)
    for k, v in annot:
        d[k].append(v.strip('"'))

    return {k: "; ".join(v) for k, v in d.items()}


def parse_gff_annot(annot):
    annot = annot.rstrip("\n").rstrip(";").split(";")
    annot = [x.strip().split("=", 1) for x in annot if x]
    # assume gff is unique
    return {k: v for k, v in annot}


def read_gtf(gtf_file, is_gff=False):
    if is_gff:
        parse_annot = parse_gff_annot
    else:
        parse_annot = parse_gtf_annot
    gene_dict = defaultdict(lambda: defaultdict(lambda: TranscriptExon()))

    # i = 0
    with open(gtf_file, "r") as f:
        for line in track(f, description="Parsing GTF..."):
            # i += 1
            # if i % 200_000 == 0:
            #     logging.info(f"Processed {i:,} lines")
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            if len(line) < 9 or line[2] != "exon":
                continue
            d = parse_annot(line[8])
            if "gene_id" in d and "transcript_id" in d:
                gene_id = d["gene_id"]
                transcript_id = d["transcript_id"]
                exon_id = d["exon_number"] if "exon_number" in d else d["exon"]
            elif "Parent" in d and "ID" in d:
                gene_id = d["Parent"]
                transcript_id, exon_id = d["ID"].rsplit(".", 1)
                exon_id = exon_id.removeprefix("exon")
            else:
                continue
            # if exon id is digit, convert to interger
            if exon_id.isdigit():
                exon_id = int(exon_id)
            if "tag" in d:
                tags = d["tag"].split("; ")
                if "MANE_Select" in tags:
                    priority = (0, 0)
                elif "Ensembl_canonical" in tags:
                    priority = (0, 1)
                elif "basic" in tags:
                    priority = (0, 2)
                else:
                    priority = (10, 0)
            elif (
                "transcript_support_level" in d
                and (sl := d.get("transcript_support_level", "").split()[0]).isdigit()
            ):
                priority = (1, int(sl))
            else:
                priority = (10, 0)
            gene_dict[gene_id][transcript_id].priority = priority
            gene_dict[gene_id][transcript_id].add_span(
                gene_id,
                d.get("gene_name", None),
                line[0],
                line[6],
                exon_id,
                int(line[3]),
                int(line[4]),
            )
    return gene_dict


def rank_exons(tx_id, exon_info):
    # level 0, MANE_Select, highest priority
    # level 1, transcript_support_level, smaller level will have higher priority
    # level 10 (default), others, lowest priority
    # patched:
    # level 2, tx_id with ".1", ".2", .. suffix, smaller number will have higher priority
    # level 3, tx_id with "-01", "-02", .. suffix, smaller number will have higher priority
    # level 4, tx_len, longer transcript will have higher priority
    if exon_info.priority[0] < 10:
        return exon_info.priority
    # some plant sample such as Arabidopsis will have ".1", ".2", ... tag in the end of tx_id
    if "." in tx_id and tx_id.split(".")[-1].isdigit():
        return (2, int(tx_id.split(".")[-1]))
    # some plant sample such as rice will have "-01" tag in the end of tx_id
    if "-" in tx_id and tx_id.split("-")[-1].isdigit():
        return (3, int(tx_id.split("-")[-1]))
    # start from 100 to 100,000, longer transcript will have higher priority
    tx_len = exon_info.calc_len()
    return (4, 100_100 - tx_len if tx_len < 100_000 else 100_001)


def parse_file(gtf_file, fasta_file, output_file, seq_file):
    gene_dict = read_gtf(
        gtf_file, is_gff=gtf_file.endswith("gff") or gtf_file.endswith("gff3")
    )
    fasta = Fasta(fasta_file, read_ahead=100_000)
    # i = 0
    with open(output_file, "w") as f1, open(seq_file, "w") as f2:
        f1.write("gene_id\ttranscript_id\tgene_name\tchrom\tstrand\tspans\n")
        for g, v in track(gene_dict.items(), description="Fetching sequences..."):
            # i += 1
            # if i % 5_000 == 0:
            #     logging.info(f"Processed {i:,} genes")
            vs = sorted(v.items(), key=lambda x: rank_exons(x[0], x[1]))
            e, v2 = vs[0]
            f1.write(f"{g}\t{e}\t{v2}\n")
            f2.write(f">{g}\n{v2.get_seq(fasta)}\n")


if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument("-g", "--gtf-file", help="GTF file", required=True)
    argparser.add_argument("-f", "--fasta-file", help="Fasta file", required=True)
    argparser.add_argument("-o", "--output-file", help="Output file", required=True)
    argparser.add_argument(
        "-s", "--seq-file", help="Output sequence file", required=True
    )

    args = argparser.parse_args()
    parse_file(args.gtf_file, args.fasta_file, args.output_file, args.seq_file)
