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

## Two-color mapping (`coralsnake.mapping`)

```python
from coralsnake.mapping import map_file, score_to_mapq

map_file(r1_file, r2_file, ref_files, output_files, unmap_file, report_file,
         forward_library, max_mismatches, threads, min_alignment_length,
         min_mapping_ratio, max_a2g_ratio, max_c2t_ratio, index_dir,
         index_only, batch_size, orientation_filter)
```

## Annotation (`coralsnake.annot`)

```python
from coralsnake.annot import run_annot

run_annot(input_file, output_file, annot_file, cols="1,2,3",
          keep_na=True, collapse_annot=False, add_count=False, skip_header=False)
```

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

## Variant analysis (`coralsnake.effect`, `.motif`, `.coordinate`)

The fused variant subcommands are importable functions.

```python
from coralsnake.effect import run_effect
from coralsnake.motif import run_motif
from coralsnake.coordinate import run_coordinate
```

Shared types/constants live in `coralsnake.effect`:

```python
from coralsnake.effect import Annot, Site, IUPAC, CODON_TABLE, reverse_base
```

- `Site` / `Annot` — input/output row dataclasses. `Annot`'s field order
  defines the output column order (identical to the standalone `variant`).
- `IUPAC` / `CODON_TABLE` / `reverse_base` — ambiguity codes, the standard
  genetic code, and strand helpers.
- `effect.build_transcript_index(gtf_file)` / `effect.site2mut(...)` — lower
  level building blocks for the classifier.

## I/O conventions

- Genomic coordinates are **0-based, half-open** internally.
- GTF is converted at the boundary (1-based closed → 0-based).
- `polars.DataFrame` is used for tabular data end to end.
