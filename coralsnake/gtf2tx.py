#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-05-01 18:23


import hashlib
from collections import defaultdict
from functools import lru_cache

import pysam
import rich.progress

from .utils import Span, Transcript, get_logger

LOGGER = get_logger(__name__)


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

    gene_dict = defaultdict(lambda: defaultdict(lambda: Transcript()))

    with rich.progress.open(gtf_file, "r", description="Parsing GTF...") as f:
        for line in f:
            if line.startswith("#"):
                continue

            line = line.strip().split("\t")
            if len(line) < 9:
                continue

            feature_type = line[2]
            if feature_type not in ["exon", "start_codon", "stop_codon"]:
                continue
            d = parse_annot(line[8])
            if (
                "gene_id" in d
                and "transcript_id" in d
                and ("exon_number" in d or "exon" in d)
            ):
                gene_id = d["gene_id"]
                transcript_id = d["transcript_id"]
                exon_id = d["exon_number"] if "exon_number" in d else d["exon"]
                if "gene_name" not in d and "gene" in d:
                    d["gene_name"] = d["gene"]
            elif "Parent" in d and "ID" in d:
                gene_id = d["Parent"]
                # eg: NbL00g00020.1.exon.1 (last 1 is the exon number)
                # eg: CsaV3_5G032770.1.exon1 (last 1 is the exon number)
                # eg: exon-XM_018814795.2-1 (last 1 is the exon number)
                # Thus, we can rsplit by ".exon." or ".exon" or "-"
                if ".exon." in d["ID"]:
                    transcript_id, exon_id = d["ID"].rsplit(".exon.", 1)
                elif ".exon" in d["ID"]:
                    transcript_id, exon_id = d["ID"].rsplit(".exon", 1)
                elif "-" in d["ID"]:
                    transcript_id, exon_id = d["ID"].rsplit("-", 1)
                    transcript_id = transcript_id.removeprefix("exon-").removeprefix(
                        "exon"
                    )
                else:
                    continue
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

            tx = gene_dict[gene_id][transcript_id]
            # do not need to check priority is set or not
            tx.priority = priority
            for k, v in zip(
                [
                    "gene_id",
                    "transcript_id",
                    "chrom",
                    "strand",
                    "gene_name",
                    "transcript_biotype",
                ],
                [
                    gene_id,
                    transcript_id,
                    line[0],
                    line[6],
                    d.get("gene_name", None),
                    d.get("transcript_biotype", None),
                ],
            ):
                if v is None:
                    continue
                if getattr(tx, k):
                    if v != getattr(tx, k):
                        raise ValueError(f"{k} mismatch")
                else:
                    setattr(tx, k, v)

            if feature_type == "exon":
                tx.add_exon(exon_id, Span(int(line[3]) - 1, int(line[4])))
            elif feature_type == "start_codon":
                tx.start_codon = Span(int(line[3]) - 1, int(line[4]))
            elif feature_type == "stop_codon":
                tx.stop_codon = Span(int(line[3]) - 1, int(line[4]))
            else:
                continue
    return gene_dict


def rank_transcript(tx_id, tx_info):
    # level 0, MANE_Select, highest priority
    # level 1, transcript_support_level, smaller level will have higher priority
    # level 10 (default), others, lowest priority
    # patched:
    # level 2, tx_id with ".1", ".2", .. suffix, smaller number will have higher priority
    # level 3, tx_id with "-01", "-02", .. suffix, smaller number will have higher priority
    # level 4, tx_len, longer transcript will have higher priority
    if tx_info.priority[0] < 10:
        return tx_info.priority
    # some plant sample such as Arabidopsis will have ".1", ".2", ... tag in the end of tx_id
    if "." in tx_id and tx_id.split(".")[-1].isdigit():
        return (2, int(tx_id.split(".")[-1]))
    # some plant sample such as rice will have "-01" tag in the end of tx_id
    if "-" in tx_id and tx_id.split("-")[-1].isdigit():
        return (3, int(tx_id.split("-")[-1]))
    # start from 100 to 100,000, longer transcript will have higher priority
    tx_len = tx_info.length
    return (4, 100_100 - tx_len if tx_len < 100_000 else 100_001)


# Define valid characters
first_char_valid_chars = set(
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$%&+./:;?@^_|~-"
)
rest_char_valid_chars = set(
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$%&*+./:;=?@^_|~-"
)


@lru_cache(maxsize=1000)
def sanitize_sequence_name(name):
    if not name:
        return name

    first_char = name[0]
    if first_char not in first_char_valid_chars:
        sanitized_first_char = "_"
    else:
        sanitized_first_char = first_char

    sanitized_rest_chars = "".join(
        c if c in rest_char_valid_chars else "_" for c in name[1:]
    )

    sanitized_name = sanitized_first_char + sanitized_rest_chars
    if sanitized_name != name:
        hash_id = hashlib.md5(name.encode()).hexdigest()[:8]
        sanitized_name = f"{sanitized_name}_{hash_id}"

    return sanitized_name


def top_transcript(gene_dict):
    for g, v in gene_dict.items():
        t, v2 = sorted(v.items(), key=lambda x: rank_transcript(x[0], x[1]))[0]
        yield v2


def parse_file(
    gtf_file,
    fasta_file,
    output_file,
    seq_file=None,
    sanitize=False,
    with_codon=False,
    with_genename=False,
    with_biotype=False,
    filter_biotype=None,
    seq_upper=False,
    line_length=0,
):
    gene_dict = read_gtf(
        gtf_file, is_gff=gtf_file.endswith("gff") or gtf_file.endswith("gff3")
    )
    if seq_file is not None:
        fasta = pysam.FastaFile(fasta_file)
        seq_writer = open(seq_file, "w")
    tsv_writer = open(output_file, "w")
    header = "gene_id\ttranscript_id\tchrom\tstrand\tspans"
    if with_codon:
        header += "\tstart_codon\tstop_codon"
    if with_genename:
        header += "\tgene_name"
    tsv_writer.write(header + "\n")
    for tx in rich.progress.track(
        top_transcript(gene_dict), total=len(gene_dict), description="Writing..."
    ):
        if filter_biotype and tx.transcript_biotype != filter_biotype:
            continue
        if sanitize:
            tx.gene_id = sanitize_sequence_name(tx.gene_id)
        if seq_file is not None:
            seq_writer.write(
                f">{tx.gene_id}\n{tx.get_seq(fasta, upper=seq_upper, wrap = line_length)}\n"
            )
        tsv_writer.write(
            tx.to_tsv(
                with_codon=with_codon,
                with_genename=with_genename,
                with_biotype=with_biotype,
            )
        )
        tsv_writer.write("\n")
    if seq_file is not None:
        seq_writer.close()
    tsv_writer.close()


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
