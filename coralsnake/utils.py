#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-25 14:21

import logging
import textwrap
from collections import defaultdict
from functools import lru_cache

import pysam

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@lru_cache(maxsize=10000)
def _python_reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


try:
    from functools import wraps

    from . import mappy as mp

    @lru_cache(maxsize=10000)
    @wraps(mp.revcomp)
    def _c_reverse_complement(seq: str) -> str:
        return mp.revcomp(seq)

    reverse_complement = _c_reverse_complement

except ImportError:
    reverse_complement = _python_reverse_complement


class Span:
    def __init__(self, start: int, end: int):
        # 0-based
        self.start = start
        self.end = end

    def __iter__(self):
        yield self.start
        yield self.end

    def __repr__(self) -> str:
        return f"(start={self.start}, end={self.end})"


class Transcript:
    def __init__(
        self,
        gene_id: str | None = None,
        transcript_id: str | None = None,
        chrom: str = "",
        strand: str = "",
        exons: dict[str | int, Span] | None = None,
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
        self._seq = None
        # extra feature
        self.gene_name: str | None = gene_name
        self.transcript_biotype: str | None = None
        self.start_codon: Span | None = None
        self.stop_codon: Span | None = None
        self.priority: tuple[int, int] = (10, 0)

    def add_exon(self, exon_id: str | int, exon: Span) -> None:
        self.exons[exon_id] = exon
        self._exons_forwards = None
        self._cum_exon_lens = None
        self._seq = None

    def sort_exons(self) -> None:
        self.exons = dict(
            sorted(
                self.exons.items(), key=lambda x: x[1].start, reverse=self.strand == "-"
            ),
        )
        self._exons_forwards = None
        self._cum_exon_lens = None
        self._seq = None

    def to_tsv(self, with_codon=False, with_genename=False, with_biotype=False) -> str:
        # convert into 1-based
        line = []
        for key in ["gene_id", "transcript_id", "chrom", "strand"]:
            v = getattr(self, key)
            if v is not None:
                line.append(v)
            else:
                line.append("")
        line += [",".join([f"{v.start + 1}-{v.end}" for v in self.exons.values()])]
        if with_codon:
            # use the start postion of the span
            line += [
                str(self.start_codon.start + 1) if self.start_codon else "",
                str(self.stop_codon.start + 1) if self.stop_codon else "",
            ]
        if with_genename:
            line.append(self.gene_name if self.gene_name is not None else "")
        if with_biotype:
            line.append(self.transcript_biotype if self.transcript_biotype else "")
        return "\t".join(line)

    def get_genome_spans(self) -> str:
        # 1-based
        return ",".join([f"{v.start + 1}-{v.end}" for v in self.exons.values()])

    def get_gene_spans(self) -> str:
        # 1-based
        gene_spans = []
        exon_start = 1
        for exon_end in self.cum_exon_lens:
            gene_spans.append(f"{exon_start}-{exon_end}")
            exon_start = exon_end + 1
        return ",".join(gene_spans)

    def get_seq(self, fasta: pysam.FastaFile, sort=True, upper=True, wrap=0):
        if self._seq is None:
            if sort:
                self.sort_exons()
            seq = ""
            for _, v in self.exons.items():
                e = fasta.fetch(self.chrom, v.start, v.end)
                if self.strand == "-":
                    e = reverse_complement(e)
                seq += e
            self._seq = seq

        seq_out = self._seq
        if upper:
            seq_out = seq_out.upper()
        if wrap > 0:
            seq_out = textwrap.fill(seq_out, wrap)
        return seq_out

    @property
    def exons_forwards(self) -> list[Span]:
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
            header = f.readline().strip("\n").split("\t")
            gene_id_idx = header.index("gene_id")
            transcript_id_idx = header.index("transcript_id")
            chrom_idx = header.index("chrom")
            strand_idx = header.index("strand")
            spans_idx = header.index("spans")
            if "gene_name" in header:
                gene_name_idx = header.index("gene_name")
            else:
                gene_name_idx = None
        else:
            (
                gene_id_idx,
                transcript_id_idx,
                chrom_idx,
                strand_idx,
                spans_idx,
                gene_name_idx,
            ) = 0, 1, 2, 3, 4, None
        for line in f:
            fields = line.strip("\n").split("\t")
            if len(fields) < 5:
                continue
            gene_id = fields[gene_id_idx]
            transcript_id = fields[transcript_id_idx]
            chrom = fields[chrom_idx]
            strand = fields[strand_idx]
            spans = fields[spans_idx]
            gene_name = fields[gene_name_idx] if gene_name_idx is not None else None
            exons = {
                idx: Span(int(start) - 1, int(end))
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
