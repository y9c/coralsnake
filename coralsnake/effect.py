#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Variant-effect annotation, fused into coralsnake from the standalone
# `variant` package (`variant effect`). pyensembl + varcode (the old deps)
# are replaced with a pure-Python classifier over coralsnake's GTF machinery
# (`metagene.load_gtf`). Naming and output column order match the original.

from dataclasses import dataclass, field

from .gtf import load_gtf
from .utils import get_logger, reverse_complement

LOGGER = get_logger(__name__)


# ---------------------------------------------------------------------------
# Shared constants / dataclasses (originally `variant` package, `variant.effect`)
# ---------------------------------------------------------------------------
# IUPAC ambiguity codes.
IUPAC = {
    "A": ["A"], "C": ["C"], "G": ["G"], "U": ["T"], "T": ["T"],
    "M": ["A", "C"], "R": ["A", "G"], "W": ["A", "T"], "S": ["C", "G"],
    "Y": ["C", "T"], "K": ["G", "T"], "V": ["A", "C", "G"],
    "H": ["A", "C", "T"], "D": ["A", "G", "T"], "B": ["C", "G", "T"],
    "N": ["G", "A", "T", "C"], ".": ["G", "A", "T", "C"], "-": ["G", "A", "T", "C"],
}

# Complement map (also handles IUPAC codes).
COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "M": "K", "R": "Y", "W": "W", "S": "S", "Y": "R", "K": "M",
    "V": "B", "H": "D", "D": "H", "B": "V", "N": "N", ".": ".", "-": "-",
}

# Standard genetic code (RNA codons -> one-letter amino acid).
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
    "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Map a generic varcode-style effect onto a coarser category.
FEATURE_MAPPER = {
    "AlternateStartCodon": "StartCodon", "StartLoss": "StartCodon",
    "StopLoss": "StopCodon", "ComplexSubstitution": "CDS",
    "Deletion": "CDS", "ExonLoss": "CDS", "FrameShiftTruncation": "CDS",
    "FrameShift": "CDS", "Insertion": "CDS", "PrematureStop": "CDS",
    "Substitution": "CDS", "Silent": "CDS",
    "ExonicSpliceSite": "SpliceSite", "IntronicSpliceSite": "SpliceSite",
    "SpliceAcceptor": "SpliceSite", "SpliceDonor": "SpliceSite",
}


def expand_base(base):
    """Expand an IUPAC code into the list of concrete bases."""
    return IUPAC.get(base.upper(), [base.upper()])


def reverse_base(base):
    """Reverse-complement a base string (e.g. a motif)."""
    try:
        return "".join(COMPLEMENT[b] for b in base)[::-1]
    except KeyError:
        return base[::-1]


@dataclass
class Site:
    """One input variant row."""
    chrom: str = "."
    pos: int = -1
    strand: str = "."
    ref: str = "-"
    alt: str = "N"
    extra: dict = field(default_factory=dict)

    def to_list(self):
        return [self.chrom, self.pos, self.strand, self.ref, self.alt]


@dataclass
class Annot:
    """The output tuple of a single variant-effect call.

    Field order defines the output column order (kept identical to the
    standalone `variant` package).
    """
    mut_type: str | None = None
    gene_type: str | None = None
    gene_name: str | None = None
    gene_pos: int | None = None
    transcript_name: str | None = None
    transcript_pos: int | None = None
    transcript_motif: str | None = None
    transcript_strand: str = "."
    coding_pos: int | None = None
    codon_ref: str | None = None
    aa_pos: int | None = None
    aa_ref: str | None = None
    distance2splice: int | None = None

    def __str__(self):
        return "\t".join(str(x) for x in vars(self).values())

    def get_values(self, as_string=False):
        values = list(vars(self).values())
        if as_string:
            return list(map(str, values))
        return values

    def get_names(self):
        return list(vars(self).keys())

    def rename_effect(self, rename_or_not=True):
        if rename_or_not:
            self.mut_type = FEATURE_MAPPER.get(str(self.mut_type), self.mut_type)
        return self


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


def _distance_to_splice(transcript, g_pos):
    """Closest distance to a splice site for an intronic position."""
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
    return "".join(CODON_TABLE.get(rna[i : i + 3], "X") for i in range(0, len(rna) - 2, 3))


def _tx_sequence(fasta, chrom, transcript):
    """Assemble the spliced transcript sequence from a genome FASTA."""
    seq = ""
    for ex in transcript["exons"]:
        piece = fasta.fetch(chrom, ex["g_start"], ex["g_end"])
        if transcript["strand"] == "-":
            piece = reverse_complement(piece)
        seq += piece
    return seq.upper()


def _codon_and_aa(transcript, t_pos, coding_pos, cds):
    """Return ``(codon_ref, aa_ref)`` for a CDS position."""
    codon_start = (coding_pos // 3) * 3
    codon_ref = cds[codon_start : codon_start + 3]
    aa_seq = _translate(cds)
    aa_idx = coding_pos // 3
    aa_ref = None
    if aa_idx < len(aa_seq):
        aa_ref = aa_seq[aa_idx]
    return codon_ref, aa_ref


# ---------------------------------------------------------------------------
# Effect classification per site
# ---------------------------------------------------------------------------
def _find_overlapping(transcripts_by_chrom, chrom, pos):
    return [tx for tx in transcripts_by_chrom.get(chrom, []) if _is_exonic(tx, pos)]


def _refine_cds_effect(codon_ref, site_ref, site_alt, coding_pos):
    """Map a CDS mutation onto a varcode-style effect name.

    ``coding_pos`` locates the mutated base within the coding sequence, so the
    substitution is applied at the exact position of the allele (not the first
    matching base). ``codon_ref`` is DNA (``T`` not ``U``) from the genome FASTA.
    """
    if len(site_ref) == 1 and len(site_alt) == 1:
        i = coding_pos % 3
        alt_codon = codon_ref[:i] + site_alt + codon_ref[i + 1 :]
        if alt_codon == codon_ref:
            return "Silent"
        # Same amino acid despite a different codon -> synonymous = silent.
        if CODON_TABLE.get(_to_rna(codon_ref)) == CODON_TABLE.get(_to_rna(alt_codon)):
            return "Silent"
        if "TAA" in alt_codon or "TAG" in alt_codon or "TGA" in alt_codon:
            return "PrematureStop"
        return "Substitution"
    # Indel: frameshift unless the net length change is a multiple of three.
    if (len(site_alt) - len(site_ref)) % 3 == 0:
        return "InFrameIndel"
    return "FrameShift"


def _to_rna(dna: str) -> str:
    """Convert a DNA codon (with T) to RNA (with U)."""
    return dna.replace("T", "U")


def _mut2eff(site, transcripts_by_chrom, fasta, strandness, pad):
    """Return the (top or all) effects for a single ``Site``."""
    chrom = str(site.chrom)
    pos = int(site.pos)
    overlaps = _find_overlapping(transcripts_by_chrom, chrom, pos)

    if not overlaps:
        return [Annot(mut_type="Intergenic")]

    results = []
    for tx in overlaps:
        t_pos = _transcript_pos(tx, pos)
        if t_pos is None:  # intronic position within a multi-exon transcript
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
                if site.ref and site.alt and site.ref != "-":
                    if codon_ref:
                        mtype = _refine_cds_effect(
                            codon_ref, site.ref, site.alt, coding_pos
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


# ---------------------------------------------------------------------------
# Effect severity ordering (originally `variant.effect_ordering`)
# ---------------------------------------------------------------------------
# Severe-to-benign ordering (higher index == more severe).
SEVERITY_ORDER = [
    "Intergenic",
    "Intragenic",
    "NoncodingTranscript",
    "IncompleteTranscript",
    "Intronic",
    "IntronicSpliceSite",
    "FivePrimeUTR",
    "ThreePrimeUTR",
    "Silent",
    "Substitution",
    "SpliceAcceptor",
    "SpliceDonor",
    "ExonicSpliceSite",
    "AlternateStartCodon",
    "ExonLoss",
    "InFrameIndel",
    "FrameShift",
    "FrameShiftTruncation",
    "Deletion",
    "Insertion",
    "ComplexSubstitution",
    "PrematureStop",
    "StartLoss",
    "StopLoss",
]

_RANK = {name: i for i, name in enumerate(SEVERITY_ORDER)}

_PRIORITY_TYPES = [
    "rRNA",
    "rRNA_pseudogene",
    "Mt_rRNA",
    "tRNA",
    "Mt_tRNA",
    "snoRNA",
    "snRNA",
    "scaRNA",
    "scRNA",
    "vault_RNA",
    "miRNA",
]


def get_top_effect(effects, pU_mode=False):
    """Return the most severe effect from a list of :class:`Annot` objects.

    ``pU_mode`` biases towards rRNA/tRNA/snoRNA-bearing effects first, mirroring
    the original behaviour.
    """
    if not effects:
        return None
    if pU_mode:
        for eff in effects:
            if getattr(eff, "gene_type", None) in _PRIORITY_TYPES:
                return eff

    def rank(eff):
        return _RANK.get(getattr(eff, "mut_type", None), 0)

    return max(effects, key=rank)


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
        top = get_top_effect(annots, pU_mode=pU_mode)
        annots = [top] if top is not None else annots
    if rename_effect:
        for a in annots:
            a.rename_effect(True)
    return annots


# ---------------------------------------------------------------------------
# CLI / library entry point
# ---------------------------------------------------------------------------
# Columns the effect engine understands, in output order.
_EFFECT_COLUMNS = ["chrom", "pos", "strand", "ref", "alt"]


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
    """Annotate every input variant with its predicted effect."""
    from xopen import xopen

    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    columns_index_mapper = dict(zip(_EFFECT_COLUMNS, columns_index))

    if reference_gtf is None:
        raise ValueError("`effect` requires --reference-gtf (or a built-in reference).")

    transcripts_by_chrom = build_transcript_index(reference_gtf)

    import pysam

    fasta_path = reference_transcript[0] if reference_transcript else None
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None

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
                body_lines = raw_lines  # every line is a variant (no header)

            output_handle.write(
                col_sep.join(input_header + Annot().get_names()) + "\n"
            )

            for raw in body_lines:
                if not raw.strip():
                    continue
                input_cols = raw.rstrip("\n").split(col_sep)
                site = _site_from_cols(input_cols, columns_index_mapper, strandness)
                annot_list = _mut2eff(site, transcripts_by_chrom, fasta, strandness, npad)
                for annot in annot_list:
                    output_handle.write(
                        col_sep.join(input_cols + annot.get_values(as_string=True)) + "\n"
                    )
    finally:
        if fasta is not None:
            fasta.close()


def _site_from_cols(input_cols, columns_index_mapper, strandness):
    """Build a :class:`Site` from parsed input columns, applying strand flip."""
    site = Site()
    for name, i in columns_index_mapper.items():
        if i < len(input_cols):
            setattr(site, name, input_cols[i])
    if site.ref and site.ref != "-":
        site.ref = site.ref.upper()
    if site.alt and site.alt != "N":
        site.alt = site.alt.upper()
    if strandness and site.strand == "-":
        # ref/alt are on the transcript's plus strand; flip for '-' strand sites.
        site.ref = reverse_base(site.ref)
        site.alt = reverse_base(site.alt)
    return site
