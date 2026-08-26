#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Shared constants + dataclasses for the migrated `variant` subcommands.
#
# FUSION NOTE: site/annotation dataclasses and IUPAC tables are migrated from
# the standalone `variant` package; logging/fasta helpers are reused from
# coralsnake.utils instead of being re-implemented.

from dataclasses import dataclass, field


# IUPAC ambiguity codes (migrated from variant.effect)
IUPAC = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "U": ["T"],
    "T": ["T"],
    "M": ["A", "C"],
    "R": ["A", "G"],
    "W": ["A", "T"],
    "S": ["C", "G"],
    "Y": ["C", "T"],
    "K": ["G", "T"],
    "V": ["A", "C", "G"],
    "H": ["A", "C", "T"],
    "D": ["A", "G", "T"],
    "B": ["C", "G", "T"],
    "N": ["G", "A", "T", "C"],
    ".": ["G", "A", "T", "C"],
    "-": ["G", "A", "T", "C"],
}

# Complement map (migrated from variant.effect)
COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "N": "N",
    ".": ".",
    "-": "-",
}

# Standard genetic code (RNA codons -> amino acid one-letter).
CODON_TABLE = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "STOP",
    "UAG": "STOP",
    "UGU": "C",
    "UGC": "C",
    "UGA": "STOP",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
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
    """The output tuple of a single variant effect call.

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


FEATURE_MAPPER = {
    "AlternateStartCodon": "StartCodon",
    "StartLoss": "StartCodon",
    "StopLoss": "StopCodon",
    "ComplexSubstitution": "CDS",
    "Deletion": "CDS",
    "ExonLoss": "CDS",
    "FrameShiftTruncation": "CDS",
    "FrameShift": "CDS",
    "Insertion": "CDS",
    "PrematureStop": "CDS",
    "Substitution": "CDS",
    "Silent": "CDS",
    "ExonicSpliceSite": "SpliceSite",
    "IntronicSpliceSite": "SpliceSite",
    "SpliceAcceptor": "SpliceSite",
    "SpliceDonor": "SpliceSite",
}


__all__ = [
    "Site",
    "Annot",
    "IUPAC",
    "COMPLEMENT",
    "CODON_TABLE",
    "FEATURE_MAPPER",
    "expand_base",
    "reverse_base",
    "motif",
    "coordinate",
    "effect",
    "effect_ordering",
]
