#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-25 14:21

import logging
from functools import lru_cache

import numpy as np

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class Exon:
    def __init__(self, start: int, end: int):
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        return f"Exon({self.start=}, {self.end=})"


class Transcript:
    def __init__(
        self,
        gene_id: str,
        transcript_id: str,
        gene_name: str,
        chrom: str,
        strand: str,
        spans: str,
    ):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.gene_name = gene_name
        self.chrom = chrom
        self.strand = strand
        self.exons = self._parse_exons(spans)
        self.cum_exon_lens = self._calculate_cum_exon_lens()
        self.length = self.cum_exon_lens[-1]

    def _parse_exons(self, spans: str) -> list[Exon]:
        # gff and gtf are 1-based, convert to 0-based
        exons = [
            Exon(int(start) - 1, int(end))
            for span in spans.split(",")
            for start, end in [span.split("-")]
        ]
        return exons if self.strand == "+" else list(reversed(exons))

    def _calculate_cum_exon_lens(self) -> np.ndarray:
        lengths = [exon.end - exon.start for exon in self.exons]
        return np.cumsum(lengths)

    def __repr__(self) -> str:
        return f"Transcript({self.gene_id=}, {self.transcript_id=}, {self.gene_name=}, {self.chrom=}, {self.strand=}, {self.exons=}, {self.length=})"


def load_annotation(annotation_file: str) -> dict[str, Transcript]:
    annot = {}
    with open(annotation_file, "r") as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split("\t")
            transcript = Transcript(*fields[:6])
            # annot[transcript.transcript_id] = transcript
            annot[transcript.gene_id] = transcript
    return annot


@lru_cache(maxsize=10000)
def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]
