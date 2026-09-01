---
layout: default
title: Coralsnake
nav_order: 1
---

# Coralsnake

![logo](./coralsnake_DNA.png)

> Exon-aware RNA analysis: splicing-aware liftover, metagene profiling and annotation

Coralsnake is an exon-aware RNA analysis pipeline — it assembles exons into
transcript references (`prepare`), **splices and joins** reads between
transcript and genome coordinates (`liftover`),
runs exon-aware metagene profiling, and annotates sites to genes/transcripts —
and bundles **sequence-logo** plotting.

Nucleotide-conversion (two-color / three-color) *mapping* is not part of
coralsnake any more; it lives in the dedicated
[`prismalign`](https://github.com/y9c/prismalign) package (pluggable backends:
bwamem, minimap2/mappy, pure-Python).

## Overview

How the subcommands fit together:

```text
           coralsnake — exon-aware RNA pipeline
      (read alignment is external: bwa / prismalign / ...)

  reads ── external mapper ──► aligned BAM
                                    │
                                    ▼
       ┌─────────────────────────────────────────────────┐
       │  prepare    GTF/GFF + genome.fa ► transcript.fa  │
       │             (exon-spliced transcript reference)  │
       │                                                  │
       │  liftover   one command, both directions:        │
       │             -d t2g   tx.bam ► genome.bam         │
       │             -d g2t   genome.bam ► tx.bam         │
       └─────────────────────────────────────────────────┘
                                    │
                   sites.tsv (chrom,pos,strand[,ref,alt])
                                    │
  ┌────────────┬────────────┬────────────┬────────────┬────────────┐
  ▼            ▼            ▼            ▼            ▼
  ┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐
  │  annotate  │ │   motif    │ │ coordinate │ │  metagene  │ │   group    │
  └────────────┘ └────────────┘ └────────────┘ └────────────┘ └────────────┘
   gene/transc    motif seq     chrom names    metagene      gene clusters
   + region +     (strand-      (UCSC ⇄        profile       + consensus
   effect         aware)        Ensembl)       (+ plot)

       ┌────────────┐
       │    logo    │   motifs ► sequence-logo image
       └────────────┘
```

## Installation

```bash
pip install coralsnake
```

Visualization commands (metagene plot, sequence logo) need the optional
`plot` extra:

```bash
pip install "coralsnake[plot]"
```

## Commands

| Command    | Description                                                       |
| ---------- | ----------------------------------------------------------------- |
| `prepare`  | Extract primary transcript from a GTF/GFF file.                   |
| `liftover` | Bidirectional BAM liftover; `--direction t2g` (default) → genome, `g2t` → transcript. |
| `annotate` | Unified site/variant annotation (region + gene/transcript/effect).|
| `group`    | Group genes and build a consensus sequence.                       |
| `metagene` | Metagene profiling across 5'UTR/CDS/3'UTR.                        |
| `logo`     | Plot a DNA/RNA sequence logo.                                     |
| `motif`    | Fetch a genomic motif around variant sites.                       |
| `coordinate` | Map chromosome names between coordinate systems.                |

## Documentation

- [Architecture & Design](architecture.md) — how the package is organized.
- [CLI Reference](cli.md) — full usage of every subcommand.
- [Python API](api.md) — import the functions directly.
- [Metagene](metagene.md) — metagene profiling details.
- [Logo](logo.md) — sequence-logo plotting details.
- [Variant](variant.md) — genomic variant analysis subcommand group.
