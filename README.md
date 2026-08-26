# Coralsnake

[![Pypi Releases](https://img.shields.io/pypi/v/coralsnake.svg)](https://pypi.python.org/pypi/coralsnake)
[![Downloads](https://pepy.tech/badge/coralsnake)](https://pepy.tech/project/coralsnake)

<p align="center">
    <picture>
        <img alt="coralsnake logo" src="https://coralsnake.yech.science/coralsnake_DNA.png" style="width: 50%">
    </picture>
</p>

Coralsnake is a transcriptome mapping toolkit. In addition to the
two-color mapping workflow (`prepare`, `map`, `liftover`, `annot`, `group`), it
now bundles the full **metagene** profiling analysis and a **sequence-logo**
plotter as first-class subcommands.

## Installation

```bash
pip install coralsnake
```

Optional support for the visualization commands (metagene profile plot and
sequence logo) requires the lightweight `plot` extra, which only pulls in
matplotlib when you need it:

```bash
pip install "coralsnake[plot]"
```

## Commands

| Command      | Description                                                        |
| ------------ | ------------------------------------------------------------------ |
| `prepare`    | Extract primary transcript from a GTF/GFF file.                    |
| `map`        | Map reads to a reference genome using BWA-MEM (two-color aware).   |
| `liftover`   | Remap transcriptome-aligned reads back to genome coordinates.      |
| `annot`      | Annotate a TSV of genomic sites with transcript positions.         |
| `group`      | Group genes and build a consensus sequence.                        |
| `metagene`   | Metagene profiling: distribution of sites across 5'UTR/CDS/3'UTR.  |
| `logo`       | Plot a DNA/RNA sequence logo (requires `coralsnake[plot]`).        |

## Metagene subcommand

`coralsnake metagene` is a full migration of the `metagene` package, built on
the high-performance `polars` + `ruranges` stack. It computes the distribution
of genomic sites relative to gene regions (5'UTR, CDS, 3'UTR) and can emit
binned statistics and a publication-ready profile plot.

```bash
# Using a built-in reference (GRCh38) or a custom GTF:
coralsnake metagene -i sites.tsv.gz -r GRCh38 -H -m 1,2,3 -w 5 \
                    -o output.tsv -s scores.tsv -p plot.png

coralsnake metagene -i sites.bed -g custom.gtf.gz -m 1,2,3 -w 5 \
                    -o output.tsv -s scores.tsv -p plot.png
```

List or download the built-in references:

```bash
coralsnake metagene --list
coralsnake metagene --download GRCh38
```

### Python API

The metagene functions are also importable directly from coralsnake:

```python
from coralsnake.metagene import (
    load_sites, load_reference, load_gtf,
    map_to_transcripts, normalize_positions, plot_profile,
)

sites = load_sites("sites.tsv.gz", with_header=True, meta_col_index=[0, 1, 2])
ref = load_reference("GRCh38")   # or load_gtf("custom.gtf.gz")
annotated = map_to_transcripts(sites, ref)
gene_bins, gene_stats, gene_splits = normalize_positions(
    annotated, split_strategy="median", bin_number=100
)
plot_profile(gene_bins, gene_splits, "metagene_plot.png")
```

## Logo subcommand

`coralsnake logo` plots a DNA/RNA sequence logo from a set of motif sequences.
The scoring engine is pure numpy; the renderer needs matplotlib (`plot` extra).

```bash
coralsnake logo -m ACGT -m ACGG -m CCGT -o logo.png
# or with per-motif weights from a file (seq\tcount)
coralsnake logo -i motifs.tsv -o logo.svg
```

```python
from coralsnake.logo import Mlogo

m = Mlogo(motifs=["ACGT", "ACGG", "CCGT"], to2bit=True)
m.plot(ax)  # requires matplotlib (plot extra)
```
