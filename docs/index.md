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
        coralsnake — exon-aware RNA analysis
   (read alignment is external: bwa / prismalign / ...)

   GTF/GFF + genome.fa ── prepare ──► transcript.fa
                                           │  (align reads against it)
   reads ── external mapper ───────────────┘
                                           ▼
                          liftover   tx.bam ⇄ genome.bam
                                           │
                                           ▼
            annotate · metagene · motif · coordinate · group
                                           │
                                           ▼
                                         logo
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
