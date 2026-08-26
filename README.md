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
| `annotate`   | Unified site/variant annotation (region + gene/transcript/effect). |
| `annot`      | Site labeling with transcript positions (legacy; see `annotate`).  |
| `effect`     | Variant effect classification (legacy; see `annotate`).            |
| `group`      | Group genes and build a consensus sequence.                        |
| `metagene`   | Metagene profiling: distribution of sites across 5'UTR/CDS/3'UTR.  |
| `logo`       | Plot a DNA/RNA sequence logo (requires `coralsnake[plot]`).        |
| `variant`    | Variant utilities (`motif`, `coordinate`).                         |

> **`annotate` is the merged successor of `annot` + `effect`.** One command, one
> GTF-based engine, one output schema: it labels a site with its
> gene/transcript/position and region, and (with a genome FASTA + ref/alt)
> classifies the variant effect (codon/AA). `annot` and `effect` are retained
> for backward compatibility.

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

The metagene functions are also importable directly from the flat modules:

```python
from coralsnake.io import load_sites, load_reference
from coralsnake.gtf import load_gtf
from coralsnake.annotation import map_to_transcripts, normalize_positions
from coralsnake.map_to_local import map_to_local
from coralsnake.plotting import plot_profile

sites = load_sites("sites.tsv.gz", with_header=True, meta_col_index=[0, 1, 2])
ref = load_reference("GRCh38")   # or load_gtf("custom.gtf.gz")
annotated = map_to_transcripts(sites, ref)
gene_bins, gene_stats, gene_splits = normalize_positions(
    annotated, split_strategy="median", bin_number=100
)
plot_profile(gene_bins, gene_splits, "metagene_plot.png")

# Map global coordinates to local transcript coordinates (strand-aware):
local = map_to_local(sites, ref, ref_id_col="transcript_id")
```

## Performance

The metagene core is built on the vectorized `polars` + `ruranges` stack, and
uses Rust-backed `ruranges` primitives instead of slow per-group Python applies:

- `map_to_transcripts` picks the best transcript per gene with a vectorized
  sort + `group_by().first()` (was `group_by().map_groups()` python apply) —
  **~20× faster** on realistic inputs.
- `map_to_local` uses `ruranges.numpy.group_cumsum` for strand-aware cumulative
  transcript offsets (was a hand-rolled `map_groups` apply) — **~7× faster**.
- `Mlogo` (sequence logo) builds its score matrix with vectorized `numpy`
  (`bincount` + codepoint lookup) — **~1.5× faster** and fixes a `0·log2(0)`
  NaN edge case.

## Logo subcommand

`coralsnake logo` plots a DNA/RNA sequence logo from a set of motif sequences.
The scoring engine is pure numpy; the renderer needs matplotlib (`plot` extra).

```bash
coralsnake logo -m ACGT -m ACGG -m CCGT -o logo.png
# or with per-motif weights from a file (seq\tcount)
coralsnake logo -i motifs.tsv -o logo.svg
```

```python
from coralsnake import Mlogo

m = Mlogo(motifs=["ACGT", "ACGG", "CCGT"], to2bit=True)
m.plot(ax)  # requires matplotlib (plot extra)
```

## Variant analysis

The `motif`, `coordinate` and `effect` commands are a migration of the
standalone `variant` package — fused into the top-level CLI with the old
`pyfaidx` / `urllib3` / `pyensembl`+`varcode` dependencies removed and
coralsnake's `pysam` + `ruranges` stack used instead. Naming and output format
are unchanged.

```bash
# Motif fetch (strand-aware, padded with N)
coralsnake motif -i sites.tsv -o motifs.tsv -f genome.fa -n 2,3 -w

# Chromosome-name mapping (UCSC ↔ Ensembl)
coralsnake coordinate -i sites.tsv -o mapped.tsv -M U2E

# Variant effect annotation (pure Python classifier on coralsnake GTF)
coralsnake effect -i variants.tsv -o effects.tsv \
                  --reference-gtf annotation.gtf \
                  --reference-transcript transcripts.fa -s -a
```

```python
from coralsnake.effect import Annot, Site, reverse_base
from coralsnake.motif import get_motif
from coralsnake.coordinate import run_coordinate
from coralsnake.effect import run_effect
```

## Documentation

- [Architecture & Design](DESIGN.md) — package layout and design decisions.
- Full docs site: <https://coralsnake.yech.science/> (see `docs/`).
