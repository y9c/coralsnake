---
layout: default
title: Python API
nav_order: 4
---

# Python API

All commands are also importable as first-class Python functions, so you can
build pipelines without shelling out.

## Core utilities (`coralsnake.utils`)

```python
from coralsnake.utils import Span, Transcript, reverse_complement, load_faidx
```

- `Span(start, end)` — a genomic span.
- `Transcript(gene_id, transcript_id, chrom, strand, exons)` — ordered exons
  with `add_exon`, `sort_exons`, `cum_exon_lens`, `length`, `get_genome_spans`,
  `get_gene_spans`, `get_seq`, `to_tsv`.
- `reverse_complement(seq)` — reverse-complement a nucleotide string.
- `load_faidx(path)` / `load_annotation(path)` — I/O helpers.

## Mapping (`coralsnake.mapping` removed)

`coralsnake.mapping` was removed in 0.0.222. Mapping is out of scope:

- **plain alignment** → `from bwamem import BwaAligner` (`bwamem` package)
- **nucleotide-conversion (two-color / three-color / 3-nt)** →
  `from prismalign import NColorMapper` (`prismalign` package, pluggable
  backends: bwamem / minimap2-mappy / pure-Python)

Feed the resulting BAM into `prepare` + `liftover` (both `--direction` options) as usual.

## Annotation (`coralsnake.annotate`)

The unified site/variant annotation tool is importable as a function:

```python
from coralsnake.annotate import run_annotate, Annotation

run_annotate(input_file, output_file, reference_gtf,
             reference_transcript="genome.fa", strandness=True)
```

`run_annotate(...)` labels a bare site (chrom,pos,strand) with its
gene/transcript/position + region and — with a genome FASTA + ref/alt — the
full variant effect (codon/AA + mut_type). `Annotation` is the output-row
dataclass whose `COLUMNS` order defines the output column order. A
precomputed `prepare` table can be passed as `annotation_table` for fast
table mode.

## Gene grouping (`coralsnake.genegroup`)

```python
from coralsnake.genegroup import group_genes, run_msa, consensus_sequence
```

## Metagene profiling (flat modules)

The full metagene pipeline is importable from the flat modules.

```python
from coralsnake.io import load_sites, load_reference  # TSV/BED/CSV → polars DataFrame
from coralsnake.gtf import prepare_exon_ref, load_gtf  # GTF → exon/codon reference (cached)
from coralsnake.annotation import (  # annotate sites → best transcript per gene
    map_to_transcripts,
    normalize_positions,       # bin positions + 5'UTR/CDS/3'UTR splits
    show_summary_stats,
)
from coralsnake.overlap import annotate_with_features, calculate_bin_statistics
from coralsnake.map_to_local import map_to_local  # global → local (spliced) coordinates
from coralsnake.plotting import plot_profile      # needs coralsnake[plot]
```

Example:

```python
import polars as pl
from coralsnake.io import load_sites
from coralsnake.gtf import load_gtf
from coralsnake.annotation import map_to_transcripts, normalize_positions

ref = load_gtf("annotation.gtf.gz")
sites = load_sites("sites.tsv.gz", with_header=True, meta_col_index=[0, 1, 2])
annotated = map_to_transcripts(sites, ref)
gene_bins, gene_stats, gene_splits = normalize_positions(
    annotated, split_strategy="median", bin_number=100
)
```

## Sequence logo (`coralsnake.logo`)

```python
from coralsnake import Mlogo

logo = Mlogo(motifs=["ACGT", "ACGG", "CCGT"], to2bit=True)
logo.plot(ax)   # requires coralsnake[plot]
```

`Mlogo` computes the score matrix with pure numpy; plotting needs matplotlib.

## Variant analysis (`coralsnake.annotate`, `.motif`, `.coordinate`)

The variant subcommands are importable functions.

```python
from coralsnake.annotate import run_annotate, Annotation  # unified site/variant annotation
from coralsnake.motif import run_motif
from coralsnake.coordinate import run_coordinate
```

- `run_annotate(...)` — label a site with gene/transcript/position + region and
  (with a genome FASTA + ref/alt) the full variant effect, in one fixed output
  schema. `Annotation`'s field order defines the output column order.
- `get_motif` / `run_motif` — strand-aware genomic motif fetch.
- `run_coordinate` — chromosome-name mapping (UCSC ↔ Ensembl).

## I/O conventions

- Genomic coordinates are **0-based, half-open** internally.
- GTF is converted at the boundary (1-based closed → 0-based).
- `polars.DataFrame` is used for tabular data end to end.
