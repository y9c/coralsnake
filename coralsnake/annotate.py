#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Unified site / variant annotation engine.
#
# This merges the logic of the two former tools:
#   * ``coralsnake annot``  – label a (chrom,pos,strand) site with the
#     containing gene / transcript and transcript-relative position.
#   * ``coralsnake effect`` – classify a variant's functional effect
#     (region, and - when a genome FASTA + ref/alt are supplied - codon,
#     amino-acid and splice context).
#
# Both are the same operation at different depths, so they share a single
# transcript index (built from a GTF once) and a single per-site classifier.
# The output schema is fixed: the first six columns are always produced, and
# the deeper `mut_*` / coding columns are filled only when applicable.

from dataclasses import dataclass

from .effect import (
    Site,
    build_transcript_index,
    _classify_exonic,
    _codon_and_aa,
    _distance_to_splice,
    _refine_cds_effect,
    _transcript_pos,
    _tx_sequence,
)
from .utils import get_logger

LOGGER = get_logger(__name__)

# Order in which a region/effect should be picked when several overlapping
# transcripts annotate the same site (benign -> severe, higher wins).
_REGION_RANK = {
    "Intergenic": 0,
    "Intronic": 1,
    "NoncodingTranscript": 2,
    "IncompleteTranscript": 3,
    "FivePrimeUTR": 4,
    "ThreePrimeUTR": 5,
    "IntronicSpliceSite": 6,
    "Silent": 7,
    "CDS": 8,
    "Substitution": 9,
    "InFrameIndel": 10,
    "FrameShift": 11,
    "PrematureStop": 12,
    "SpliceSite": 13,
}

_REGION_NAMES = (
    "Intergenic",
    "Intronic",
    "FivePrimeUTR",
    "CDS",
    "ThreePrimeUTR",
    "NoncodingTranscript",
)


@dataclass
class Annotation:
    """One annotation row (a site x one overlapping transcript).

    ``COLUMNS`` fixes the output column order and is identical whether the
    input was a bare site or a full variant - deeper columns are left as
    ``None`` when the data needed to fill them was not supplied.
    """

    gene_id: str | None = None
    transcript_id: str | None = None
    transcript_pos: int | None = None
    region: str | None = None
    gene_pos: int | None = None
    transcript_strand: str = "."
    mut_type: str | None = None
    transcript_motif: str | None = None
    coding_pos: int | None = None
    codon_ref: str | None = None
    aa_pos: int | None = None
    aa_ref: str | None = None
    distance2splice: int | None = None

    COLUMNS = [
        "gene_id",
        "transcript_id",
        "transcript_pos",
        "region",
        "gene_pos",
        "transcript_strand",
        "mut_type",
        "transcript_motif",
        "coding_pos",
        "codon_ref",
        "aa_pos",
        "aa_ref",
        "distance2splice",
    ]

    def to_list(self):
        return [getattr(self, c) for c in self.COLUMNS]

    @classmethod
    def header(cls):
        return list(cls.COLUMNS)


def _find_overlapping(transcripts_by_chrom, chrom, pos):
    """Return transcripts whose exonic span contains ``pos``."""
    return [
        tx for tx in transcripts_by_chrom.get(chrom, []) if _is_exonic(tx, pos)
    ]


def _is_exonic(transcript, g_pos):
    return any(ex["g_start"] <= g_pos < ex["g_end"] for ex in transcript["exons"])


def _annotate_site(
    site,
    transcripts_by_chrom,
    fasta=None,
    pad=10,
    strandness=True,
):
    """Classify one site into a list of :class:`Annotation` (one per overlap).

    ``fasta`` is optional; without it (or without ref/alt), the coding-effect
    columns (motif / codon / amino-acid) are left as ``None`` but the region
    and gene/transcript position are always computed.
    """
    chrom = str(site.chrom)
    pos = int(site.pos)
    overlaps = _find_overlapping(transcripts_by_chrom, chrom, pos)

    if not overlaps:
        return [Annotation(region="Intergenic", mut_type="Intergenic")]

    results = []
    for tx in overlaps:
        tid = tx.get("transcript_id")
        t_pos = _transcript_pos(tx, pos)

        if t_pos is None:
            # Inside the gene body but not in any exon -> intronic.
            results.append(
                Annotation(
                    gene_id=tx["gene_id"],
                    transcript_id=tid,
                    region="Intronic",
                    mut_type="Intronic",
                    transcript_strand=tx["strand"],
                    distance2splice=_distance_to_splice(tx, pos),
                )
            )
            continue

        region = _classify_exonic(tx, t_pos)
        gene_pos = t_pos + (tx["start_codon_pos"] or 0)
        transcript_motif = None
        coding_pos = codon_ref = aa_pos = aa_ref = None
        tx_seq = None

        if fasta is not None and tx["chrom"] in fasta.references:
            tx_seq = _tx_sequence(fasta, tx["chrom"], tx)
            s5 = tx_seq[max(t_pos - pad, 0) : t_pos].rjust(pad, "N")
            s0 = tx_seq[t_pos] if t_pos < len(tx_seq) else "N"
            s3 = tx_seq[t_pos + 1 : t_pos + 1 + pad].ljust(pad, "N")
            transcript_motif = s5 + s0 + s3

        mut_type = region
        if region == "CDS":
            start = tx["start_codon_pos"]
            stop = tx["stop_codon_pos"]
            if start is not None and fasta is not None and tx_seq is not None:
                cds = tx_seq[start : (stop + 1) if stop is not None else None]
                coding_pos = t_pos - start
                codon_ref, aa_ref = _codon_and_aa(tx, t_pos, coding_pos, cds)
                aa_pos = coding_pos // 3 + 1
                if site.ref and site.alt and site.ref != "-":
                    if codon_ref:
                        mut_type = _refine_cds_effect(
                            codon_ref, site.ref, site.alt, coding_pos
                        )
            else:
                # Coding region but no FASTA / incomplete annotation on hand.
                mut_type = "IncompleteTranscript" if fasta is not None else "CDS"

        results.append(
            Annotation(
                gene_id=tx["gene_id"],
                transcript_id=tid,
                transcript_pos=t_pos,
                region=region,
                gene_pos=gene_pos,
                transcript_strand=tx["strand"],
                mut_type=mut_type,
                transcript_motif=transcript_motif,
                coding_pos=coding_pos,
                codon_ref=codon_ref,
                aa_pos=aa_pos,
                aa_ref=aa_ref,
                distance2splice=_distance_to_splice(tx, pos),
            )
        )

    if not results:
        results = [Annotation(region="Intergenic", mut_type="Intergenic")]
    return results


def _pick_top(annotations):
    """Return the single most severe annotation (for --top / default mode)."""
    if not annotations:
        return None
    return max(
        annotations, key=lambda a: _REGION_RANK.get(a.mut_type or a.region, 0)
    )


# ---------------------------------------------------------------------------
# Library / CLI entry point
# ---------------------------------------------------------------------------
_ANNOTATE_COLUMNS = ["chrom", "pos", "strand", "ref", "alt"]


def run_annotate(
    input_file,
    output_file,
    reference_gtf,
    reference_transcript=None,
    npad=10,
    strandness=True,
    all_effects=False,
    with_header=False,
    columns="1,2,3,4,5",
):
    """Annotate every input row (site or variant) with the unified schema."""
    from xopen import xopen

    if reference_gtf is None:
        raise ValueError("`annotate` requires --reference-gtf (or a built-in reference).")

    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    columns_index_mapper = dict(zip(_ANNOTATE_COLUMNS, columns_index))

    transcripts_by_chrom = build_transcript_index(reference_gtf)

    import pysam

    fasta = (
        pysam.FastaFile(reference_transcript[0]) if reference_transcript else None
    )

    try:
        with xopen(output_file, "wt") as output_handle:
            with xopen(input_file, "rt") as input_handle:
                raw_lines = input_handle.readlines()

            if with_header:
                header_line = raw_lines[0].rstrip("\n")
                body_lines = raw_lines[1:]
                input_header = header_line.split(col_sep)
            else:
                first = raw_lines[0].rstrip("\n").split(col_sep)
                input_header = ["."] * len(first)
                for n, i in columns_index_mapper.items():
                    if i < len(input_header):
                        input_header[i] = n
                body_lines = raw_lines

            output_handle.write(
                col_sep.join(
                    input_header
                    + [
                        "gene_id",
                        "transcript_id",
                        "transcript_pos",
                        "region",
                        "gene_pos",
                        "transcript_strand",
                        "mut_type",
                        "transcript_motif",
                        "coding_pos",
                        "codon_ref",
                        "aa_pos",
                        "aa_ref",
                        "distance2splice",
                    ]
                )
                + "\n"
            )

            for raw in body_lines:
                if not raw.strip():
                    continue
                input_cols = raw.rstrip("\n").split(col_sep)
                site = _site_from_cols(input_cols, columns_index_mapper, strandness)
                annotations = _annotate_site(
                    site, transcripts_by_chrom, fasta, npad, strandness
                )
                if not all_effects:
                    top = _pick_top(annotations)
                    annotations = [top] if top is not None else annotations
                for ann in annotations:
                    out_cols = ["" if v is None else str(v) for v in ann.to_list()]
                    output_handle.write(
                        col_sep.join(input_cols + out_cols) + "\n"
                    )
    finally:
        if fasta is not None:
            fasta.close()


def _site_from_cols(input_cols, columns_index_mapper, strandness):
    """Build a :class:`Site` from parsed input columns, applying strand flip."""
    from .effect import reverse_base

    site = Site()
    for name, i in columns_index_mapper.items():
        if i < len(input_cols):
            value = input_cols[i]
            if name == "pos":
                try:
                    value = int(value)
                except ValueError:
                    pass
            setattr(site, name, value)
    if site.ref and site.ref != "-":
        site.ref = site.ref.upper()
    if site.alt and site.alt != "N":
        site.alt = site.alt.upper()
    if strandness and site.strand == "-":
        site.ref = reverse_base(site.ref) if site.ref else "-"
        site.alt = reverse_base(site.alt) if site.alt not in ("N", "") else site.alt
    return site


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("input_file")
    ap.add_argument("output_file")
    ap.add_argument("reference_gtf")
    ap.add_argument("--reference-transcript", "-f", default=None, nargs="+")
    ap.add_argument("--columns", "-c", default="1,2,3,4,5")
    ap.add_argument("--all-effects", "-a", action="store_true")
    ap.add_argument("--npad", "-n", type=int, default=10)
    ap.add_argument("--with-header", action="store_true")
    args = ap.parse_args()
    run_annotate(
        args.input_file,
        args.output_file,
        args.reference_gtf,
        args.reference_transcript,
        args.npad,
        True,
        args.all_effects,
        args.with_header,
        args.columns,
    )
