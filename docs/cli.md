---
layout: default
title: CLI Reference
nav_order: 2
---

# CLI Reference

```
Usage: coralsnake [OPTIONS] COMMAND [ARGS]...

  Coralsnake (transcriptome mapping utils)
```

## Commands

| Command    | Description                                                       |
| ---------- | ----------------------------------------------------------------- |
| `prepare`  | Extract primary transcript from gtf/gff file.                     |
| `map`      | Map reads to reference genome using BWA-MEM.                      |
| `liftover` | Fetch genomic motif.                                              |
| `annot`    | Annotate tsv file.                                                |
| `group`    | Group and find consensus of gene.                                 |
| `metagene` | Run metagene profiling analysis on genomic sites.                 |
| `logo`     | Plot a DNA/RNA sequence-logo (requires `coralsnake[plot]`).       |

Run `coralsnake <command> --help` for the full option list of any command.
