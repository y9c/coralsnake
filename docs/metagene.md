---
layout: default
title: Metagene
nav_order: 5
---

# Metagene

`coralsnake metagene` is a full migration of the `metagene` package, built on
the high-performance `polars` + `ruranges` stack. It computes the distribution
of genomic sites relative to gene regions (5'UTR, CDS, 3'UTR) and can emit
binned statistics and a publication-ready profile plot.

## Usage

```bash
# Using a built-in reference (GRCh38):
coralsnake metagene -i sites.tsv.gz -r GRCh38 -H -m 1,2,3 -w 5 \
                    -o output.tsv -s scores.tsv -p plot.png

# Using a custom GTF:
coralsnake metagene -i sites.bed -g custom.gtf.gz -m 1,2,3 -w 5 \
                    -o output.tsv -s scores.tsv -p plot.png
```

List or download the built-in references (use the `reference` command; the
`metagene --list` / `--download` flags still work but are deprecated):

```bash
coralsnake reference list
coralsnake reference download GRCh38
```

## Key options

| Option              | Description                                                     |
| ------------------- | --------------------------------------------------------------- |
| `-i, --input`       | Input file path (BED, TSV, CSV, ...).                           |
| `-o, --output`      | Output annotated TSV/CSV.                                       |
| `-s, --output-score`| Output binned score statistics.                                 |
| `-p, --output-figure`| Output metagene profile plot (requires `coralsnake[plot]`).     |
| `-r, --reference`   | Built-in reference (e.g. `GRCh38`, `GRCm39`).                   |
| `-g, --gtf`         | GTF/GFF file for a custom reference.                            |
| `-m, --meta-columns`| 1-based column indices for Chromosome,Start[,End],Strand.       |
| `-w, --weight-columns`| 1-based column indices for weight/score values.               |
| `-b, --bins`        | Number of bins (default 100).                                   |
| `-H, --with-header` | Input has a header line.                                        |

## Python API

```python
from coralsnake.io import load_sites, load_reference
from coralsnake.gtf import load_gtf
from coralsnake.annotation import map_to_transcripts, normalize_positions
from coralsnake.map_to_local import map_to_local
from coralsnake.plotting import plot_profile

sites = load_sites("sites.tsv.gz", with_header=True, meta_col_index=[0, 1, 2])
ref = load_reference("GRCh38")          # or load_gtf("custom.gtf.gz")
annotated = map_to_transcripts(sites, ref)
gene_bins, gene_stats, gene_splits = normalize_positions(
    annotated, split_strategy="median", bin_number=100
)
plot_profile(gene_bins, gene_splits, "metagene_plot.png")

local = map_to_local(sites, ref, ref_id_col="transcript_id")
```

## Performance

The metagene core uses Rust-backed `ruranges` primitives rather than slow
per-group Python applies:

- `map_to_transcripts`: vectorized best-transcript selection (`~20×` faster).
- `map_to_local`: `ruranges.numpy.group_cumsum` for cumsum offsets (`~7×` faster).
