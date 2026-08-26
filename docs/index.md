---
layout: default
title: Coralsnake
nav_order: 1
---

# Coralsnake

![logo](./coralsnake_DNA.png)

> Transcriptome mapping in two colors

Coralsnake is a transcriptome mapping toolkit with a two-color mapping
workflow and bundled **metagene** profiling and **sequence-logo** plotting.

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
| `map`      | Map reads to a reference genome using BWA-MEM (two-color aware).  |
| `liftover` | Remap transcriptome-aligned reads back to genome coordinates.     |
| `annot`    | Annotate a TSV of genomic sites with transcript positions.        |
| `group`    | Group genes and build a consensus sequence.                       |
| `metagene` | Metagene profiling across 5'UTR/CDS/3'UTR.                        |
| `logo`     | Plot a DNA/RNA sequence logo.                                     |
| `motif`    | Fetch a genomic motif around variant sites.                       |
| `coordinate` | Map chromosome names between coordinate systems.                |
| `effect`   | Annotate genomic variant effects.                                 |

## Documentation

- [Architecture & Design](architecture.md) — how the package is organized.
- [CLI Reference](cli.md) — full usage of every subcommand.
- [Python API](api.md) — import the functions directly.
- [Metagene](metagene.md) — metagene profiling details.
- [Logo](logo.md) — sequence-logo plotting details.
- [Variant](variant.md) — genomic variant analysis subcommand group.
