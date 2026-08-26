#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Migrated from the standalone `variant` package (`variant effect`).
#
# pyensembl + varcode (the old, buggy deps) are replaced with a pure-Python
# effect classifier built on coralsnake's GTF machinery (`metagene.load_gtf`),
# which already provides transcript-relative exon offsets and start/stop
# codon positions. Naming and output column order match the standalone
# `variant` package.

import sys

from ..metagene import load_gtf
from ..utils import get_logger, reverse_complement
from . import Annot, CODON_TABLE, Site
from . import effect_ordering

LOGGER = get_logger(__name__)


# ---------------------------------------------------------------------------
# Transcript index helpers
# ---------------------------------------------------------------------------
def _transcripts_from_gtf(reference):
    """Build a per-transcript structure from the polars exon reference.

    ``reference`` is the output of :func:`coralsnake.metagene.load_gtf`.
    Returns ``{transcript_id: {...}}`` with 0-based transcript-relative
    cumulative exon offsets, strand, and start/stop codon positions.
    """
    import polars as pl

    tx = {}
    for (chrom, strand), grp in reference.group_by(
        ["Chromosome", "Strand"], maintain_order=True
    ):
        grp = grp.sort("Start", descending=(strand == "-"))
        lengths = pl.col("End") - pl.col("Start")
        grp = grp.with_columns(
            lengths.cum_sum().over("transcript_id").alias("_cum_ce")
        ).with_columns(
            (lengths.cum_sum().over("transcript_id") - lengths).alias("_cum_cs")
        )
        acc = {}
        for r in grp.iter_rows(named=True):
            tid = r["transcript_id"]
            if tid not in acc:
                acc[tid] = {
                    "strand": strand,
                    "chrom": chrom,
                    "gene_id": r["gene_id"],
                    "transcript_length": r["transcript_length"],
                    "start_codon_pos": r["start_codon_pos"],
                    "stop_codon_pos": r["stop_codon_pos"],
                    "exons": [],
                }
            acc[tid]["exons"].append(
                {
                    "g_start": int(r["Start"]),
                    "g_end": int(r["End"]),
                    "t_start": int(r["_cum_cs"]),
                    "t_end": int(r["_cum_ce"]),
                }
            )
        tx.update(acc)
    return tx


def build_transcript_index(gtf_file: str):
    """Load a GTF into ``{chromosome: [transcript, ...]}`` for the classifier."""
    reference = load_gtf(gtf_file)
    txs = _transcripts_from_gtf(reference)
    by_chrom = {}
    for _, tx in txs.items():
        by_chrom.setdefault(tx["chrom"], []).append(tx)
    return by_chrom


# ---------------------------------------------------------------------------
# Position helpers (strand-aware, 0-based)
# ---------------------------------------------------------------------------
def _transcript_pos(transcript, g_pos):
    """Map a genomic position to a 0-based transcript-relative position.

    Returns ``None`` when the position is intronic / outside all exons.
    """
    for ex in transcript["exons"]:
        if ex["g_start"] <= g_pos < ex["g_end"]:
            if transcript["strand"] == "+":
                return g_pos - ex["g_start"] + ex["t_start"]
            return ex["g_end"] - 1 - g_pos + ex["t_start"]
    return None


def _is_exonic(transcript, g_pos):
    return any(ex["g_start"] <= g_pos < ex["g_end"] for ex in transcript["exons"])


def _classify_exonic(transcript, t_pos):
    """Classify an exonic position by effect type."""
    start = transcript["start_codon_pos"]
    stop = transcript["stop_codon_pos"]
    if start is None:
        return "NoncodingTranscript"
    if t_pos < start:
        return "FivePrimeUTR"
    if stop is not None and t_pos > stop:
        return "ThreePrimeUTR"
    return "CDS"


def _classify_intronic(transcript, g_pos):
    """Classify an intronic position (between exons)."""
    exons = transcript["exons"]
    for i in range(len(exons) - 1):
        left, right = exons[i], exons[i + 1]
        if left["g_end"] <= g_pos < right["g_start"]:
            dist_donor = g_pos - left["g_end"]
            dist_acceptor = right["g_start"] - 1 - g_pos
            if transcript["strand"] == "-":
                dist_donor, dist_acceptor = dist_acceptor, dist_donor
            if dist_donor < 2:
                return "SpliceDonor"
            if dist_acceptor < 2:
                return "SpliceAcceptor"
            return "Intronic"
    return "Intronic"


def _distance_to_splice(transcript, g_pos):
    exons = transcript["exons"]
    if len(exons) < 2:
        return None
    d2s = []
    for i, ex in enumerate(exons):
        if i == 0:
            if transcript["strand"] == "+":
                d2s.append(g_pos - ex["g_end"])
            else:
                d2s.append(ex["g_start"] - g_pos)
        elif i == len(exons) - 1:
            if transcript["strand"] == "+":
                d2s.append(g_pos - ex["g_start"])
            else:
                d2s.append(ex["g_end"] - g_pos)
        else:
            if transcript["strand"] == "+":
                d2s.append(g_pos - ex["g_start"])
                d2s.append(g_pos - ex["g_end"])
            else:
                d2s.append(ex["g_start"] - g_pos)
                d2s.append(ex["g_end"] - g_pos)
    if not d2s:
        return None
    return int(sorted(d2s, key=abs)[0])


# ---------------------------------------------------------------------------
# Sequence + codon helpers
# ---------------------------------------------------------------------------
def _translate(cds: str) -> str:
    rna = cds.replace("T", "U")
    return "".join(
        CODON_TABLE.get(rna[i : i + 3], "X") for i in range(0, len(rna) - 2, 3)
    )


def _tx_sequence(fasta, chrom, transcript):
    seq = ""
    for ex in transcript["exons"]:
        piece = fasta.fetch(chrom, ex["g_start"], ex["g_end"])
        if transcript["strand"] == "-":
            piece = reverse_complement(piece)
        seq += piece
    return seq.upper()


def _codon_and_aa(transcript, t_pos, coding_pos, cds):
    """Return (codon_ref, aa_ref) for a CDS position."""
    codon_start = (coding_pos // 3) * 3
    codon_ref = cds[codon_start : codon_start + 3]
    aa_ref = None
    aa_seq = _translate(cds)
    aa_idx = coding_pos // 3
    if aa_idx < len(aa_seq):
        aa_ref = aa_seq[aa_idx]
    return codon_ref, aa_ref


# ---------------------------------------------------------------------------
# Effect classification per site
# ---------------------------------------------------------------------------
def _find_overlapping(transcripts_by_chrom, chrom, pos):
    return [tx for tx in transcripts_by_chrom.get(chrom, []) if _is_exonic(tx, pos)]


def _mut2eff(site, transcripts_by_chrom, fasta, strandness, pad):
    chrom = str(site.chrom)
    pos = int(site.pos)
    overlaps = _find_overlapping(transcripts_by_chrom, chrom, pos)

    if not overlaps:
        # Intergenic: only report the gene nearest the site.
        return [Annot(mut_type="Intergenic", gene_name=None)]

    results = []
    for tx in overlaps:
        t_pos = _transcript_pos(tx, pos)
        is_exonic = t_pos is not None
        if not is_exonic:
            # Shouldn't happen given _is_exonic, but guard anyway.
            continue

        mtype = _classify_exonic(tx, t_pos)
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

        if mtype == "CDS":
            start = tx["start_codon_pos"]
            stop = tx["stop_codon_pos"]
            if start is not None and fasta is not None and tx_seq is not None:
                cds = tx_seq[start : (stop + 1) if stop is not None else None]
                coding_pos = t_pos - start
                codon_ref, aa_ref = _codon_and_aa(tx, t_pos, coding_pos, cds)
                aa_pos = coding_pos // 3 + 1
                # Refine generic CDS into a specific substitution name.
                if site.ref and site.alt and site.ref != "-":
                    if codon_ref:
                        mtype = _refine_cds_effect(
                            mtype, codon_ref, site.ref, site.alt, aa_ref
                        )
            else:
                mtype = "IncompleteTranscript"

        results.append(
            Annot(
                mut_type=mtype,
                gene_type=None,
                gene_name=tx["gene_id"],
                gene_pos=gene_pos,
                transcript_name=tx.get("name", tx["gene_id"]),
                transcript_pos=t_pos,
                transcript_motif=transcript_motif,
                transcript_strand=tx["strand"],
                coding_pos=coding_pos,
                codon_ref=codon_ref,
                aa_pos=aa_pos,
                aa_ref=aa_ref,
                distance2splice=_distance_to_splice(tx, pos),
            )
        )

    if not results:
        return [Annot(mut_type="Intergenic")]
    return results


def _refine_cds_effect(mtype, codon_ref, site_ref, site_alt, aa_ref):
    """Map a CDS mutation onto a varcode-style effect name.

    This is intentionally conservative: it distinguishes substitutions that
    are silent, missense, or introduce a stop codon. Insertions/deletions are
    reported as frameshifts when the change is not a multiple of three.

    ``codon_ref`` is DNA (``T`` not ``U``) since it comes from the genome FASTA.
    """
    if len(site_ref) == 1 and len(site_alt) == 1:
        # Single-nucleotide substitution within a codon.
        alt_codon = codon_ref.replace(site_ref, site_alt, 1)
        if alt_codon == codon_ref:
            return "Silent"
        # DNA stop codons.
        if "TAA" in alt_codon or "TAG" in alt_codon or "TGA" in alt_codon:
            return "PrematureStop"
        return "Substitution"
    # Indel: frameshift unless length change is a multiple of three.
    if (len(site_alt) - len(site_ref)) % 3 == 0:
        return "InFrameIndel"
    return "FrameShift"


def site2mut(
    site,
    transcripts_by_chrom,
    fasta=None,
    pad=10,
    strandness=True,
    all_effects=False,
    rename_effect=False,
    pU_mode=False,
):
    """Back-compat entry point mirroring the standalone ``site2mut``."""
    annots = _mut2eff(site, transcripts_by_chrom, fasta, strandness, pad)
    if not all_effects and len(annots) > 1:
        # Pick the top effect using the migrated effect ordering.
        top = effect_ordering.get_top_effect(annots, pU_mode=pU_mode)
        annots = [top] if top is not None else annots
    if rename_effect:
        for a in annots:
            a.rename_effect(True)
    return annots


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def run_effect(
    input_file,
    output_file,
    reference_gtf,
    reference_transcript,
    reference_protein,
    npad,
    strandness,
    all_effects,
    pU_mode,
    with_header,
    columns,
):
    import gzip

    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    columns_index_mapper = dict(
        zip(["chrom", "pos", "strand", "ref", "alt"], columns_index)
    )

    if reference_gtf is None:
        LOGGER.error("`variant effect` requires --reference-gtf (or a built-in).")
        sys.exit(1)

    transcripts_by_chrom = build_transcript_index(reference_gtf)

    import pysam

    fasta_path = reference_transcript[0] if reference_transcript else None
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None

    def open_in(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

    def open_out(path):
        return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

    with open_out(output_file) as output_handle:
        with open_in(input_file) as input_handle:
            raw_lines = input_handle.readlines()

        # Determine header + body lines.
        if with_header:
            header_line = raw_lines[0].rstrip("\n")
            body_lines = raw_lines[1:]
            input_header = header_line.split(col_sep)
            header = col_sep.join(input_header + Annot().get_names()) + "\n"
        else:
            first = raw_lines[0].rstrip("\n").split(col_sep)
            input_header = ["."] * len(first)
            for n, i in columns_index_mapper.items():
                if i < len(input_header):
                    input_header[i] = n
            header = col_sep.join(input_header + Annot().get_names()) + "\n"
            body_lines = raw_lines  # process every line (no header)

        output_handle.write(header)

        for raw in body_lines:
            if not raw.strip():
                continue
            input_cols = raw.rstrip("\n").split(col_sep)
            site = Site()
            for n, i in columns_index_mapper.items():
                if i < len(input_cols):
                    setattr(site, n, input_cols[i])
            if site.ref and site.ref != "-":
                site.ref = site.ref.upper()
            if site.alt and site.alt != "N":
                site.alt = site.alt.upper()
            if strandness and site.strand == "-":
                # ref/alt are on the plus strand of the transcript; flip if site strand is '-'
                from . import reverse_base

                site.ref = reverse_base(site.ref)
                site.alt = reverse_base(site.alt)
            annot_list = _mut2eff(site, transcripts_by_chrom, fasta, strandness, npad)
            for annot in annot_list:
                output_handle.write(
                    col_sep.join(input_cols + annot.get_values(as_string=True)) + "\n"
                )

    if fasta is not None:
        fasta.close()
