---
layout: default
title: Architecture & Design
nav_order: 2
---

# Coralsnake Architecture & Design

This document explains how coralsnake is organized, why it is designed the way
it is, and how the pieces fit together for anyone extending it.

> Lifecycle note: coralsnake grew from a single-purpose two-color transcriptome
> mapper into an exon-aware RNA pipeline; the mapping engine itself now lives in
> the [`prismalign`](https://github.com/y9c/prismalign) package (on top of the
> lightweight [`bwamem`](https://github.com/y9c/bwamem) binding). This page is the
> canonical "why is it built this way" reference.

## 1. Package layout

```
coralsnake/
├── cli.py              # click/rich_click command group — the only entry point
├── gtf2tx.py           # GTF/GFF → reference transcript FASTA
├── tbam2gbam.py        # transcriptome BAM → genome BAM (liftover)
├── gbam2tbam.py        # reverse direction (thin wrapper)
├── genegroup.py        # gene clustering & consensus
├── annotate.py         # unified site/variant annotation (annot + effect)
├── annot.py            # legacy table-based site annotation
├── utils.py            # core data structures (Span, Transcript) + logging/plot/io helpers
├── seqops.c            # C kernels: RC, batch conversion, score_and_tag, reverse_md
├── logo.py             # DNA/RNA sequence-logo (numpy scoring, optional mpl)
├── motif.py            # genomic motif fetch (strand-aware)
├── coordinate.py       # chrom-name mapping (UCSC↔Ensembl)
├── effect.py           # variant-effect annotation + Annot/Site, genetic code
├── annotation.py       # metagene: map_to_transcripts, normalize_positions
├── map_to_local.py     # metagene: global→local transcript coordinates
├── overlap.py          # metagene: feature overlap + bin statistics
├── io.py               # metagene: site / reference loading
├── gtf.py              # metagene: GTF → exon/codon reference (cached)
├── download.py         # metagene: builtin-reference downloader
├── config.py           # metagene: builtin reference registry
├── plotting.py         # metagene: profile plot (optional matplotlib)
└── __init__.py         # lazy top-level exports (Mlogo, __version__)
```

## 2. The command layer

`coralsnake.cli:cli` is a `rich_click` group. Each subcommand is a thin wrapper
that validates input and then delegates to a module function. The rule of thumb:

> **CLI does argument parsing; modules do the work.** Keeping the CLI thin makes
> every function importable and testable from Python.

| Command    | Module function        | Purpose                                    |
| ---------- | ---------------------- | ------------------------------------------ |
| `prepare`  | `gtf2tx.parse_file`    | Extract primary transcript from GTF/GFF.    |
| `liftover --direction t2g` | `tbam2gbam.convert_bam` | Remap transcriptome BAM to genome BAM. |
| `liftover --direction g2t` | `gbam2tbam.convert_bam` | Remap genome BAM to transcriptome BAM. |
| `annotate` | `annotate.run_annotate`| Unified site/variant annotation.             |
| `group`    | `genegroup.group_genes`| Cluster genes & build consensus sequences.  |
| `metagene` | `annotation.*`         | Metagene profiling across 5'UTR/CDS/3'UTR.  |
| `motif`    | `motif.run_motif`      | Strand-aware genomic motif fetch.           |
| `coordinate`| `coordinate.run_coordinate` | Chrom-name mapping (UCSC↔Ensembl).    |
| `logo`     | `logo.Mlogo`           | Plot a DNA/RNA sequence logo.               |

> `map` was removed in 0.0.222 — mapping lives in
> [`prismalign`](https://github.com/y9c/prismalign) (nucleotide-conversion /
> two-color) and [`bwamem`](https://github.com/y9c/bwamem) (plain BWA-MEM,
> including its `HierarchicalAligner` for priority-ordered multi-reference
> mapping).

## 3. Layered design

### 3.1 The performance core
The hot paths are all built on high-performance vectorized primitives rather
than per-row Python loops:

- **ruranges** (Rust-backed) for genomic interval operations:
  `overlaps`, `group_cumsum`, `spliced_subsequence`, `merge`, `nearest`, etc.
- **polars** (Rust-backed dataframe) for table transformations and aggregations.
- **numpy** for array math and binned statistics.

The explicit design rule is: **never iterate intervals in Python inside a hot
path**. Where a `group_by().map_groups()` (Python apply) appeared in the
original code, it was replaced with a vectorized equivalent (sort + window
functions, or a ruranges primitive).

Example — best-transcript selection in `annotation.py`:
```python
# Before (slow): group_by().map_groups(python_apply)  → ~20x slower
# After (fast):  sort by (level asc, length desc, txid asc) then group_by().first()
annot = annot.sort(["gene_id", "transcript_level", "transcript_length", "transcript_id"],
                   descending=[False, False, True, False])
best_tx = annot.group_by("gene_id", maintain_order=True).first()
annot = annot.join(best_tx, on=["gene_id", "transcript_id"])
```

### 3.2 Strand conventions
Coralsnake is strand-aware throughout. A single convention is followed
consistently:

- **Genomic coordinates are 0-based, half-open** `[start, end)` internally.
- GTF files arrive 1-based closed and are converted at the boundary.
- ruranges' `group_cumsum` uses `negative_strand=True` to mean **plus** strand
  (`strand == "+"`) — this is documented in `map_to_local._strand_aware_cumsum`.

### 3.3 Coordinate systems
| Function | Input | Output |
| -------- | ----- | ------ |
| `map_to_local` | genome coords | transcript (spliced) coords |
| `spliced_subsequence` | local slice | genome coords |
| `remap_to_genome` | transcriptome-aligned BAM | genome BAM |

## 4. The metagene profiling pipeline

The metagene profiling pipeline is a set of flat modules (`gtf.py`, `io.py`,
`annotation.py`, `map_to_local.py`, `overlap.py`, `plotting.py`) exposed as the
`coralsnake metagene` subcommand and an importable Python API. Its data flow:

```
GTF ──load_gtf──▶ exon/codon reference (polars) ──┐
sites ──load_sites──▶ site table (polars) ────────┤
                                                  ▼
                                        map_to_transcripts
                                                  │  (per-site: best transcript)
                                                  ▼
                                        normalize_positions
                                                  │  (bins + 5'UTR/CDS/3'UTR splits)
                                                  ▼
                                    gene_bins / gene_stats / gene_splits
                                                  │
                                  ┌───────────────┴───────────────┐
                                  ▼                               ▼
                               plot_profile               write_csv (score table)
                          (optional matplotlib)          (feature_type + bin stats)
```

Built-in references (GRCh38, mm39, ...) are registered in `config.py` and
downloaded on demand to the XDG cache dir by `download.py`. The whole pipeline
works offline with a user-supplied GTF.

## 5. Optional dependencies

Coralsnake keeps the base install light. matplotlib is the only heavy optional
dependency and is gated behind the `plot` extra:

```python
# coralsnake/logo.py
def _require_plotting():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError as e:
        raise ImportError("Plotting requires 'coralsnake[plot]'") from e
```

Both the metagene profile plot and the sequence logo lazily import matplotlib,
so the core analysis works without it. `coralsnake[plot]` adds it back.

## 6. Key design decisions

1. **`ruranges.numpy.overlaps`, not top-level `ruranges.overlaps`.** The modern
   ruranges API places interval kernels under `ruranges.numpy`. Standardize on
   it.
2. **Preserve byte-identical output on refactors.** Every vectorized rewrite was
   verified against the previous behavior (throughout the migration) before
   landing — see the README performance section.
3. **Centralize I/O.** `xopen` handles gzip transparently for FASTA/TSV; `pysam`
   handles BAM/FASTA; `polars` handles tabular I/O.
4. **Cache aggressively.** `load_gtf` writes a sidecar `.parquet` (and a pickle
   cache in `annot`) so repeated runs on the same reference are fast.

## 7. The `variant` commands (`motif`, `coordinate`)

The standalone [`variant`](https://github.com/y9c/variant) toolkit has been
fused into coralsnake as flat top-level modules (no `variant` subpackage), with
its old buggy dependencies removed:

| Command | Old dep (removed)   | New implementation |
| ------- | ------------------- | ------------------ |
| `motif` | pyfaidx | pysam.FastaFile |
| `coordinate` | urllib3 | stdlib urllib |

Naming and output column order are kept identical to the standalone package.
The variant-effect classifier (region + codon/AA annotation, formerly
`pyensembl` + `varcode`) is now exposed through the unified `annotate` command,
which reuses coralsnake's GTF/CDS machinery (transcript offsets, start/stop
codon positions) instead of re-implementing it — a conservative pure-Python
classifier lives in `annotate.py`.

### 7.1 Duplicate-fusion policy
Where the standalone package re-implemented things coralsnake already had, the
migration **fused** rather than duplicated:

- `reverse_complement`, `reverse_base`, `expand_base` → reuse coralsnake
  `utils.reverse_complement` / `variant` IUPAC tables.
- GTF → transcript exon/codon indexing → reuse `metagene.load_gtf`.
- gzip/streaming I/O → reuse the same pattern as `motif`/`coordinate`.
- `Site` / `Annot` dataclasses and their field order → kept verbatim so
  output formats remain byte-compatible with the old package.

