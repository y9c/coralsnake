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


def _transcript_span(transcript):
    """Genomic (0-based half-open) span from first to last exon."""
    return (
        min(e["g_start"] for e in transcript["exons"]),
        max(e["g_end"] for e in transcript["exons"]),
    )


def _find_overlapping(transcripts_by_chrom, chrom, pos):
    """Return transcripts whose exon OR gene-body span contains ``pos``.

    Including the body span lets an intronic position (inside the gene but not
    in any exon) be classified as Intronic rather than Intergenic.
    """
    hits = []
    for tx in transcripts_by_chrom.get(chrom, []):
        if not tx["exons"]:
            continue
        lo, hi = _transcript_span(tx)
        if _is_exonic(tx, pos) or lo <= pos < hi:
            hits.append(tx)
    return hits


def _is_exonic(transcript, g_pos):
    return any(ex["g_start"] <= g_pos < ex["g_end"] for ex in transcript["exons"])


def build_span_index(transcripts_by_chrom):
    """Per-chromosome sorted gene-body spans for fast interval lookup.

    Returns ``{chrom: [(body_start, body_end, transcript), ...]}`` sorted by
    body_start, so a site can be looked up with ``bisect`` instead of scanning
    every transcript on the chromosome.
    """
    index = {}
    for chrom, txs in transcripts_by_chrom.items():
        spans = []
        for tx in txs:
            if not tx["exons"]:
                continue
            lo, hi = _transcript_span(tx)
            spans.append((lo, hi, tx))
        spans.sort(key=lambda x: x[0])  # by body_start
        index[chrom] = spans
    return index


def _overlap_transcripts_fast(span_index, chrom, pos):
    """Transcripts whose gene body (or an exon) contains ``pos``, via bisect."""
    import bisect

    spans = span_index.get(chrom)
    if not spans:
        return []
    starts = [s[0] for s in spans]
    n = bisect.bisect_right(starts, pos)  # first body_start > pos
    hits = [tx for (lo, hi, tx) in spans[:n] if hi > pos]
    return hits


def _annotate_site(
    site,
    transcripts_by_chrom,
    fasta=None,
    pad=10,
    strandness=True,
    span_index=None,
    seq_cache=None,
):
    """Classify one site into a list of :class:`Annotation` (one per overlap).

    ``fasta`` is optional; without it (or without ref/alt), the coding-effect
    columns (motif / codon / amino-acid) are left as ``None`` but the region
    and gene/transcript position are always computed. ``span_index`` (from
    ``build_span_index``) avoids a full per-site scan of every transcript.
    ``seq_cache`` (optional per-run dict) lets us assemble each transcript's
    spliced sequence once instead of re-fetching from the FASTA for every site -
    a large speedup when many sites fall in the same gene.
    """
    chrom = str(site.chrom)
    pos = int(site.pos)
    if span_index is not None:
        overlaps = _overlap_transcripts_fast(span_index, chrom, pos)
    else:
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
            if seq_cache is not None:
                tx_seq = seq_cache.get(tid)
                if tx_seq is None:
                    tx_seq = _tx_sequence(fasta, tx["chrom"], tx)
                    seq_cache[tid] = tx_seq
            else:
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

# Unified output columns (see Annotation.COLUMNS) - kept in sync.
_UNIFIED_COLUMNS = [
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
    annotation_table=None,
):
    """Annotate every input row with the unified schema.

    Two input modes are supported by the *same* engine:
      * GTF mode  - ``reference_gtf`` (+ optional ``reference_transcript`` FASTA).
        Yields region + gene/transcript position + (with ref/alt + FASTA) the
        full variant effect.
      * Table mode - ``annotation_table`` (a precomputed `prepare`-style table
        with gene_id/transcript_id/chrom/strand/spans). A fast vectorised path
        yielding gene/transcript/transcript_pos (the legacy `annot` behaviour);
        region-dependent columns are left empty.
    """
    from xopen import xopen

    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    columns_index_mapper = dict(zip(_ANNOTATE_COLUMNS, columns_index))

    if annotation_table is not None:
        return _run_annotate_table(
            input_file,
            output_file,
            annotation_table,
            columns,
            with_header,
        )

    if reference_gtf is None:
        raise ValueError(
            "`annotate` requires --reference-gtf (or --annotation / --reference-transcript)."
        )

    transcripts_by_chrom = build_transcript_index(reference_gtf)
    span_index = build_span_index(transcripts_by_chrom)

    import pysam

    fasta = (
        pysam.FastaFile(reference_transcript[0]) if reference_transcript else None
    )
    # Per-run cache of assembled transcript sequences (motif/codon path).
    seq_cache = {} if fasta is not None else None

    def _emit(site, input_cols, output_handle):
        annotations = _annotate_site(
            site, transcripts_by_chrom, fasta, npad, strandness, span_index, seq_cache
        )
        if not all_effects:
            top = _pick_top(annotations)
            annotations = [top] if top is not None else annotations
        for ann in annotations:
            out_cols = ["" if v is None else str(v) for v in ann.to_list()]
            output_handle.write(col_sep.join(input_cols + out_cols) + "\n")

    def _build_header(first):
        input_header = ["."] * len(first)
        for n, i in columns_index_mapper.items():
            if i < len(input_header):
                input_header[i] = n
        return input_header

    try:
        with xopen(output_file, "wt") as output_handle, xopen(input_file, "rt") as input_handle:
            first = input_handle.readline()
            if not first:
                return  # empty input: nothing to annotate
            if with_header:
                input_header = first.rstrip("\n").split(col_sep)
            else:
                input_header = _build_header(first.rstrip("\n").split(col_sep))
            output_handle.write(col_sep.join(input_header + _UNIFIED_COLUMNS) + "\n")

            if not with_header:
                input_cols = first.rstrip("\n").split(col_sep)
                try:
                    site = _site_from_cols(input_cols, columns_index_mapper, strandness)
                    _emit(site, input_cols, output_handle)
                except (ValueError, IndexError):
                    pass  # skip a malformed first row

            for raw in input_handle:
                if not raw.strip():
                    continue
                input_cols = raw.rstrip("\n").split(col_sep)
                try:
                    site = _site_from_cols(input_cols, columns_index_mapper, strandness)
                    _emit(site, input_cols, output_handle)
                except (ValueError, IndexError):
                    # Skip a malformed row instead of aborting the whole run.
                    continue
    finally:
        if fasta is not None:
            fasta.close()


def _run_annotate_table(input_file, output_file, annot_file, columns, skip_header):
    """Fast precomputed-table mode (subsumes the legacy ``annot`` command).

    Uses the same vectorised table engine as ``coralsnake.annot`` but emits the
    unified ``annotate`` schema (gene_id/transcript_id/transcript_pos filled;
    region/effect columns empty).
    """
    from .annot import _annotate_batch, _read_sites, parse_annot_file

    tree, info = parse_annot_file(annot_file, cache=True)
    cols = [int(i) - 1 for i in str(columns).split(",")]
    from xopen import xopen

    with xopen(input_file, "rt") as fi, xopen(output_file, "wt") as fo:
        lines, chroms, positions, strands = _read_sites(fi, cols, skip_header)
        results = _annotate_batch(lines, chroms, positions, strands, tree, info)
        write_header = True
        for line, annot_list in zip(lines, results):
            if write_header:
                fo.write(line + "\tgene_id\ttranscript_id\ttranscript_pos\t"
                         "region\tgene_pos\ttranscript_strand\tmut_type\t"
                         "transcript_motif\tcoding_pos\tcodon_ref\taa_pos\taa_ref\t"
                         "distance2splice\n")
                write_header = False
            if annot_list:
                for gene_id, transcript_id, transcript_pos in annot_list:
                    fo.write(
                        f"{line}\t{gene_id}\t{transcript_id}\t{transcript_pos}"
                        f"\t\t\t\t\t\t\t\t\t\n"
                    )
            else:
                fo.write(f"{line}\t.\t.\t.\t\t\t\t\t\t\t\t\t\n")


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
                    raise ValueError(
                        f"Position column not an integer: {value!r}"
                    )
            setattr(site, name, value)
    # A usable site needs at least a chromosome and a position.
    if "chrom" not in columns_index_mapper or "pos" not in columns_index_mapper:
        raise ValueError("columns must identify 'chrom' and 'pos'")
    if columns_index_mapper["chrom"] >= len(input_cols) or columns_index_mapper["pos"] >= len(
        input_cols
    ):
        raise ValueError("row is missing required chromosome/position columns")
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
