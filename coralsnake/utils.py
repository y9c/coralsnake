#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-25 14:21

import logging
from collections import defaultdict
from functools import lru_cache

from pyfaidx import Fasta

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@lru_cache(maxsize=10000)
def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]


class Exon:
    def __init__(self, start: int, end: int):
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        return f"Exon(start={self.start}, end={self.end})"


class Transcript:
    def __init__(
        self,
        gene_id: str | None = None,
        transcript_id: str | None = None,
        chrom: str = "",
        strand: str = "",
        exons: dict[str | int, Exon] | None = None,
        gene_name: str | None = None,
    ):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.exons = {} if exons is None else exons
        # calculated feature
        self._exons_forwards = None
        self._cum_exon_lens = None
        # extra feature
        self.gene_name = gene_name
        self.priority: tuple[int, int] = (10, 0)

    def add_exon(self, exon_id: str | int, exon: Exon) -> None:
        self.exons[exon_id] = exon
        self._exons_forwards = None
        self._cum_exon_lens = None

    def sort_exons(self) -> None:
        self.exons = dict(
            sorted(
                self.exons.items(), key=lambda x: x[1].start, reverse=self.strand == "-"
            ),
        )
        self._exons_forwards = None
        self._cum_exon_lens = None

    def to_tsv(self) -> str:
        # convert into 1-based
        line = []
        spans = ",".join([f"{v.start + 1}-{v.end}" for v in self.exons.values()])
        for key in ["gene_id", "transcript_id", "gene_name", "chrom", "strand"]:
            v = getattr(self, key)
            if v is not None:
                line.append(v)
            else:
                line.append("")
        return "\t".join(line + [spans])

    def get_seq(self, fasta: Fasta, sort=True):
        if sort:
            self.sort_exons()
        seq = ""
        for _, v in self.exons.items():
            e = fasta[self.chrom][v.start : v.end]
            if self.strand == "-":
                e = e.reverse.complement
            seq += e.seq
        return seq.upper()

    @property
    def exons_forwards(self) -> list[Exon]:
        if self._exons_forwards is None:
            self._exons_forwards = list(self.exons.values())
            if self.strand == "-":
                self._exons_forwards = self._exons_forwards[::-1]
        return self._exons_forwards

    @property
    def cum_exon_lens(self) -> list[int]:
        if self._cum_exon_lens is None:
            lengths = [exon.end - exon.start for exon in self.exons_forwards]
            cum_lengths = []
            total = 0
            for length in lengths:
                total += length
                cum_lengths.append(total)
            return cum_lengths
        return self._cum_exon_lens

    @property
    def length(self) -> int:
        return (self.cum_exon_lens or [0])[-1]

    def __repr__(self) -> str:
        res = []
        for key in [
            "gene_id",
            "transcript_id",
            "gene_name",
            "chrom",
            "strand",
            # "exons_forwards",
        ]:
            res.append(f"{key}={getattr(self, key)}")
        return f"Transcript({', '.join(res)})"


def load_annotation(
    annotation_file: str, with_header: bool = True
) -> dict[str, dict[str, Transcript]]:
    annot = defaultdict(dict)
    with open(annotation_file, "r") as f:
        if with_header:
            next(f)  # Skip header
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue
            gene_id, transcript_id, gene_name, chrom, strand, spans = fields[:6]
            exons = {
                idx: Exon(int(start) - 1, int(end))
                for idx, span in enumerate(spans.split(","), 1)
                for start, end in [span.split("-")]
            }
            transcript = Transcript(
                gene_id=gene_id,
                transcript_id=transcript_id,
                chrom=chrom,
                strand=strand,
                exons=exons,
                gene_name=gene_name,
            )
            annot[gene_id][transcript_id] = transcript
    return annot


def load_faidx(faidx_file: str) -> dict[str, int]:
    faidx = {}
    with open(faidx_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            faidx[fields[0]] = int(fields[1])
    return faidx
