# Coralsnake

[![Pypi Releases](https://img.shields.io/pypi/v/coralsnake.svg)](https://pypi.python.org/pypi/coralsnake)
[![Downloads](https://pepy.tech/badge/coralsnake)](https://pepy.tech/project/coralsnake)

<p align="center">
    <picture>
        <img alt="coralsnake logo" src="https://coralsnake.yech.science/coralsnake_DNA.png" style="width: 50%">
    </picture>
</p>

Coralsnake is an **exon-aware RNA analysis pipeline**. Its core is the
exon structure of RNA: it assembles exons into transcript references
(`prepare`), and then **splices and joins** reads between transcript and
genome coordinates (`liftover`), runs exon-aware
metagene profiling, and annotates sites to genes/transcripts. It also bundles a
sequence-logo plotter as a first-class subcommand.

Nucleotide-conversion (two-color / three-color) mapping is **not** part of
coralsnake any more: it lives in the dedicated
[`prismalign`](https://github.com/y9c/prismalign) package (pluggable backends:
bwamem, minimap2/mappy, pure-Python), on top of the lightweight
[`bwamem`](https://github.com/y9c/bwamem) BWA-MEM binding.

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

Optional support for the visualization commands (metagene profile plot and
sequence logo) requires the lightweight `plot` extra, which only pulls in
matplotlib when you need it:

```bash
pip install "coralsnake[plot]"
```

## Commands

| Route | Command | Description |
| ----- | ------- | ----------- |
| t -> g | `liftover` | Remap transcriptome-aligned reads to genome coordinates (default `--direction t2g`). |
| g -> t | `liftover -d g2t` | Remap genome-aligned reads back to transcript coordinates. |

## Read mapping (both directions)

`prepare` builds a transcript reference. The BAM-conversion commands
round-trip between transcript and genome span:
- `coralsnake liftover` (default `--direction t2g`) – transcript BAM → genome
  BAM (splices reads at exon boundaries, inserts introns).
- `coralsnake liftover -d g2t` – genome BAM → transcript BAM (clips to exons,
  joins spliced reads contiguously on the transcript).

> **Mapping itself is out of scope:** align reads with `bwamem map`, with
> `prismalign map` for nucleotide-conversion (two/three-color) chemistry, or any
> other mapper, then feed the BAM into the commands above wearing a matching
> reference.

## Command reference
| Command      | Description                                                        |
| ------------ | ------------------------------------------------------------------ |
| `prepare`    | Extract primary transcript from a GTF/GFF file.                    |
| `liftover`  | Remap reads between genome/transcript coords (`--direction t2g` default, `g2t` inverts). |
| `annotate`   | Unified site/variant annotation (region + gene/transcript/effect). |
| `group`      | Group genes and build a consensus sequence.                        |
| `metagene`   | Metagene profiling: distribution of sites across 5'UTR/CDS/3'UTR.  |
| `logo`       | Plot a DNA/RNA sequence logo (requires `coralsnake[plot]`).        |
| `motif`      | Fetch a genomic motif around variant sites (strand-aware).         |
| `coordinate` | Map chromosome names between coordinate systems (UCSC↔Ensembl).    |

> **`annotate` is the single annotation tool** — one command, one schema, two
> input modes sharing one engine:
> - `--reference-gtf [--reference-transcript FASTA]` – region + gene/transcript/
>   position + (with ref/alt + FASTA) the full variant effect.
> - `--annotation <prepare-table>` – fast precomputed-table site labeling.

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

The `motif` and `coordinate` commands are a migration of the standalone
`variant` package — fused into the top-level CLI with the old `pyfaidx` /
`urllib3` dependencies removed and coralsnake's `pysam` + `ruranges` stack used
instead. Naming and output format are unchanged. Variant *effect* annotation is
covered by `annotate` (`--reference-gtf` + a genome FASTA + ref/alt columns).

```bash
# Motif fetch (strand-aware, padded with N)
coralsnake motif -i sites.tsv -o motifs.tsv -f genome.fa -n 2,3 -w

# Chromosome-name mapping (UCSC ↔ Ensembl)
coralsnake coordinate -i sites.tsv -o mapped.tsv -M U2E

# Variant effect annotation (region + codon/AA + effect)
coralsnake annotate -i variants.tsv -o effects.tsv \
                    --reference-gtf annotation.gtf \
                    --reference-transcript genome.fa -s -a
```

```python
from coralsnake.motif import get_motif
from coralsnake.coordinate import run_coordinate
from coralsnake.annotate import run_annotate
```

## Documentation

- [Architecture & Design](DESIGN.md) — package layout and design decisions.
- Full docs site: <https://coralsnake.yech.science/> (see `docs/`).
