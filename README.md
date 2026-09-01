# Coralsnake

[![Pypi Releases](https://img.shields.io/pypi/v/coralsnake.svg)](https://pypi.python.org/pypi/coralsnake)
[![Downloads](https://pepy.tech/badge/coralsnake)](https://pepy.tech/project/coralsnake)

<p align="center">
    <picture>
        <img alt="coralsnake logo" src="https://coralsnake.yech.science/coralsnake_DNA.png" style="width: 50%">
    </picture>
</p>

RNA analysis is hard because RNA is structured three ways at once: an
**abundance hierarchy** (ribosomal rRNA is orders of magnitude more abundant
than mRNA), **strands** (sense vs. antisense), and **splicing** (mRNAs are
assembled from exons, so a transcript does not line up with the genome).
Coralsnake is an **exon-aware RNA analysis pipeline** built around these
properties: it turns a GTF/GFF into spliced transcript references (`prepare`),
**splices and joins** reads between transcript and genome coordinates in both
directions (`liftover`), and runs the analyses you do on the results —
**annotate** places sites/variants on the RNA hierarchy (5'UTR / CDS / 3'UTR /
intronic / intergenic) and calls the variant effect, **metagene** profiles how
sites distribute across 5'UTR / CDS / 3'UTR, **motif** fetches the strand-aware
reference motif around each site, and **logo** renders a DNA/RNA sequence logo.

<img src="https://coralsnake.yech.science/coralsnake_overview.svg" alt="Coralsnake pipeline overview" style="width: 720px; max-width: 100%;" />

## Installation

```bash
pip install coralsnake
```

The visualization commands (metagene plot, sequence logo) need the lightweight
`plot` extra, which only pulls in matplotlib when you need it:

```bash
pip install "coralsnake[plot]"
```

## Commands

| Command      | What it does                                                        |
| ------------ | ------------------------------------------------------------------- |
| `prepare`    | Extract the spliced primary transcript reference from GTF/GFF.      |
| `liftover`   | Splice/join reads between genome/transcript BAMs (`-d t2g` default, `-d g2t` inverts). |
| `annotate`   | Unified site/variant annotation: region on the RNA hierarchy + gene/transcript + variant effect. |
| `metagene`   | Exon-aware metagene profiling across 5'UTR/CDS/3'UTR.              |
| `motif`      | Strand-aware genomic motif fetch around variant sites.              |
| `coordinate` | Map chromosome names between coordinate systems (UCSC↔Ensembl).    |
| `group`      | Group genes and build a consensus sequence.                         |
| `logo`       | Plot a DNA/RNA sequence logo (needs `coralsnake[plot]`).            |

> **`annotate` is the single annotation tool** — one command, one schema, two
> input modes: `--reference-gtf` (region + gene/transcript, and the full variant
> effect when given a genome FASTA + ref/alt) or `--annotation <table>` (fast
> precomputed-table site labeling).

## Quick example

A typical end-to-end run (read alignment is done by any external mapper, e.g.
`bwa` / `prismalign`):

```bash
# 1. Build the spliced transcript reference from the annotation
coralsnake prepare -g annotation.gtf -f genome.fa -o transcript.fa --with-txpos

# 2. Align reads to transcript.fa with an external mapper  →  tx.bam

# 3. Splice the transcript-aligned BAM back to genome coordinates
coralsnake liftover -d t2g -i tx.bam -o genome.bam -a annotation.tsv -f genome.fai

# 4. Annotate sites to genes/transcripts (RNA hierarchy + variant effect)
coralsnake annotate -i sites.tsv -o annotated.tsv \
                    --reference-gtf annotation.gtf \
                    --reference-transcript genome.fa -s -a

# 5. Exon-aware metagene profile across 5'UTR / CDS / 3'UTR
coralsnake metagene -i sites.tsv -g annotation.gtf -o profile.tsv -p profile.png
```

## Subcommands

### `prepare` — spliced transcript reference

Build the spliced transcript reference (the target that reads are aligned to)
from a GTF/GFF and a genome FASTA:

```bash
coralsnake prepare -g annotation.gtf -f genome.fa -o transcript.fa \
                   --with-codon --with-genename --filter-biotype protein_coding
```

### `liftover` — splice-aware BAM conversion

`prepare` builds the transcript reference; `liftover` round-trips a BAM between
transcript and genome coordinates, splicing reads at exon boundaries:

- `coralsnake liftover -d t2g` (default) — transcript BAM → genome BAM (splices
  reads at exon boundaries, inserts introns).
- `coralsnake liftover -d g2t` — genome BAM → transcript BAM (clips to exons,
  joins spliced reads contiguously on the transcript).

### `annotate` — exon-aware annotation

Label a site (chrom,pos,strand) with gene/transcript/position + region
(5'UTR/CDS/3'UTR/intronic/intergenic); add a genome FASTA + ref/alt to get the
full variant effect (codon/AA + mut_type). Fast precomputed-table mode via
`--annotation <table>`.

```bash
coralsnake annotate -i sites.tsv -o out.tsv --reference-gtf annotation.gtf -c 1,2,3
coralsnake annotate -i variants.tsv -o effects.tsv \
                    --reference-gtf annotation.gtf \
                    --reference-transcript genome.fa -s -a
```

### `metagene` — exon-aware metagene profiling

Built on the high-performance `polars` + `ruranges` stack. Computes the
distribution of sites relative to gene regions (5'UTR, CDS, 3'UTR) and can emit
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

**Python API** — the functions are importable from the flat modules:

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

### `motif` & `coordinate`

A migration of the standalone `variant` package — fused into the top-level CLI
with the old `pyfaidx` / `urllib3` dependencies removed and coralsnake's
`pysam` + `ruranges` stack used instead. Naming and output format are unchanged.

```bash
# Motif fetch (strand-aware, padded with N)
coralsnake motif -i sites.tsv -o motifs.tsv -f genome.fa -n 2,3 -w

# Chromosome-name mapping (UCSC ↔ Ensembl)
coralsnake coordinate -i sites.tsv -o mapped.tsv -M U2E
```

```python
from coralsnake.motif import get_motif
from coralsnake.coordinate import run_coordinate
from coralsnake.annotate import run_annotate
```

### `group` — gene clustering & consensus

Group related genes and build a consensus sequence:

```bash
coralsnake group -f genes.fa -g genes.gtf -o grouped.tsv \
                 --output-consensus consensus.fa --threads 8
```

### `logo` — sequence logo

Plots a DNA/RNA sequence logo from a set of motif sequences. The scoring engine
is pure numpy; the renderer needs matplotlib (`plot` extra).

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

## Performance

The core is built on the vectorized `polars` + `ruranges` stack, using Rust-backed
`ruranges` primitives instead of slow per-group Python applies:

- `map_to_transcripts` picks the best transcript per gene with a vectorized
  sort + `group_by().first()` (was `group_by().map_groups()` python apply) —
  **~20× faster** on realistic inputs.
- `map_to_local` uses `ruranges.numpy.group_cumsum` for strand-aware cumulative
  transcript offsets (was a hand-rolled `map_groups` apply) — **~7× faster**.
- `Mlogo` (sequence logo) builds its score matrix with vectorized `numpy`
  (`bincount` + codepoint lookup) — **~1.5× faster** and fixes a `0·log2(0)`
  NaN edge case.

## Documentation

- [Architecture & Design](DESIGN.md) — package layout and design decisions.
- Full docs site: <https://coralsnake.yech.science/> (see `docs/`).

---

> **Mapping is out of scope.** Nucleotide-conversion (two-color / three-color)
> mapping is not part of coralsnake any more: it lives in the dedicated
> [`prismalign`](https://github.com/y9c/prismalign) package (pluggable backends:
> bwamem, minimap2/mappy, pure-Python), on top of the lightweight
> [`bwamem`](https://github.com/y9c/bwamem) BWA-MEM binding.
