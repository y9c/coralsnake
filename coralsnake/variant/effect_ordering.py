#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Migrated from the standalone `variant` package (`variant.effect_ordering`).
# The original relied on varcode's EffectCollection rankings; we reimplement
# the same severity ordering over our own lightweight Annot objects.

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

    # pU mode: prefer effects whose gene/biotype is a priority RNA type.
    if pU_mode:
        for eff in effects:
            if getattr(eff, "gene_type", None) in _PRIORITY_TYPES:
                return eff

    def rank(eff):
        return _RANK.get(getattr(eff, "mut_type", None), 0)

    return max(effects, key=rank)
