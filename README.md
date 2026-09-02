# Coralsnake

[![Pypi Releases](https://img.shields.io/pypi/v/coralsnake.svg)](https://pypi.python.org/pypi/coralsnake)
[![Downloads](https://pepy.tech/badge/coralsnake)](https://pepy.tech/project/coralsnake)

<p align="center">
    <picture>
        <img alt="coralsnake logo" src="https://coralsnake.yech.science/coralsnake_DNA.png" style="width: 50%">
    </picture>
</p>

Conventional genomics tools are built around DNA — a linear, double-stranded
reference where each locus lines up 1:1 with the sequence — and often lack
explicit consideration of RNA's structural properties: an **abundance
hierarchy** (ribosomal rRNA is orders of magnitude more abundant than
mRNA), **strand orientation** (sense vs. antisense), and **splicing** (mRNAs
are assembled from exons, so a transcript does not line up with its genomic
locus). Coralsnake is an **exon-aware RNA analysis pipeline** built around
exactly these properties: it cleans and normalizes the reference inputs
(`refine`), turns a GTF/GFF into spliced transcript references (`prepare`),
**splices and joins** reads between transcript and genome coordinates in both
directions (`liftover`), and runs the analyses you do on the results —
**annotate** places sites/variants on the RNA hierarchy (5'UTR / CDS / 3'UTR /
intronic / intergenic) and calls the variant effect, **metagene** profiles how
sites distribute across 5'UTR / CDS / 3'UTR, **motif** fetches the strand-aware
reference motif around each site, and **logo** renders a DNA/RNA sequence logo.

<img src="https://coralsnake.yech.science/coralsnake_overview.svg?v=2" alt="Coralsnake pipeline overview" style="width: 720px; max-width: 100%;" />

## Installation

Requires Python ≥ 3.12.

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
| `refine`     | Clean/normalize the genome FASTA + GTF before `prepare` (rename seqnames, normalize names/biotypes, flag canonical transcripts). |
| `prepare`    | Extract the spliced primary transcript reference from GTF/GFF.      |
| `liftover`   | Splice/join reads between genome/transcript BAMs (`-d t2g` default, `-d g2t` inverts). |
| `annotate`   | Unified site/variant annotation: region on the RNA hierarchy + gene/transcript + variant effect. |
| `metagene`   | Exon-aware metagene profiling across 5'UTR/CDS/3'UTR (profile plot needs `coralsnake[plot]`). |
| `motif`      | Strand-aware genomic motif fetch around variant sites.              |
| `coordinate` | Map chromosome names between coordinate systems (UCSC↔Ensembl).    |
| `group`      | Group genes and build a consensus sequence.                         |
| `logo`       | Plot a DNA/RNA sequence logo, or export its score matrix via `--matrix` (plotting needs `coralsnake[plot]`). |

> **`annotate` is the single annotation tool** — one command, one schema, two
> input modes: `--reference-gtf` (region + gene/transcript, and the full variant
> effect when given a genome FASTA + ref/alt) or `--annotation <table>` (fast
> precomputed-table site labeling).

## Quick example

A typical end-to-end run (read alignment is done by any external mapper, e.g.
`bwa` / `prismalign`):

```bash
# 0. (Optional) Clean/normalize the FASTA + GTF first (rename seqnames,
#    normalize names/biotypes, flag the canonical transcript). When you run
#    this, feed outdir/ref.annotation.gtf / outdir/ref.genome.fasta into the
#    steps below instead of annotation.gtf / genome.fa.
coralsnake refine -f genome.fa -g annotation.gtf -o outdir -n ref

# 1. Build the spliced transcript reference: FASTA to align to (-s) and the
#    annotation table (-o) used by liftover and annotate
coralsnake prepare -g annotation.gtf -f genome.fa \
                   -s transcript.fa -o annotation.tsv --with-codon

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

### `refine` — reference cleaning (pre-`prepare`)

Clean and normalize a genome FASTA and/or GTF before building the spliced
reference: seqname renaming/filtering, gene/transcript name and type
normalization, missing gene/transcript/exon row creation, overlapping-exon
merge, canonical-transcript flagging, and coordinate checks. Codon and UTR
features are preserved, so the output stays usable for `metagene` and
`annotate`. Indexing uses pysam — no external samtools/bgzip/tabix.

The GTF side is built on the shared `coralsnake.genemodel.GeneModel` object —
the same read/serialize layer `prepare` uses — so a refined GTF is a drop-in
replacement for the input: both GENCODE-style (`*_type`) and Ensembl-style
(`*_biotype`) biotype attributes are written, and the selected canonical
transcript is tagged `Ensembl_canonical` so `prepare`'s ranking picks the
transcript `refine` selected.

```bash
coralsnake refine -f genome.fa -g annotation.gtf -o outdir -n hg38 \
                  -m chrom_map.tsv -c canonical.tsv
```

- `-f/--fasta-file` — the genome FASTA to clean: seqnames are renamed/filtered,
  headers rewritten as `>new_name old_name`, and a `.fai` + `.genome.sizes`
  index is rebuilt with pysam.
- `-g/--gtf-file` — the GTF to normalize (the GeneModel pipeline above).
- `-o/--outdir` — output directory (default `./`).
- `-n/--name` — output name prefix (default: the outdir's basename).
- `-m/--rename-mapper` — TSV mapping old seqname → new seqname.
- `-p/--seqname-pattern` — keep only seqnames matching this regex.
- `-c/--canonical-transcripts` — TSV of canonical transcript IDs (1st column).
- At least one of `-f`/`-g` is required. Pass both to keep the FASTA and the
  annotation seqnames in sync — the same `-m`/`-p` apply to both. Outputs land
  under `<outdir>/<name>.{genome.fasta,annotation.gtf}` (+ `.gz`, `.fai`,
  `.genome.sizes`, tabix index), plus `.skip.gtf` (genes that failed checks)
  and `.gene_features_summary.txt`.

### `prepare` — spliced transcript reference

Build the spliced transcript reference from a GTF/GFF and a genome FASTA.
Two outputs:

- `-s/--seq-file` — the spliced transcript **FASTA**, the target that reads
  are aligned to (requires `-f`).
- `-o/--output-file` — the annotation **table** (TSV: gene/transcript, chrom,
  strand, spliced exon spans, plus optional codon/genename/biotype/txpos
  columns), consumed by `liftover -a`, `liftover --table`, and
  `annotate --annotation`.

```bash
coralsnake prepare -g annotation.gtf -f genome.fa \
                   -s transcript.fa -o annotation.tsv \
                   --with-codon --with-genename --filter-biotype protein_coding
```

### `liftover` — splice-aware BAM conversion

`prepare` builds the transcript reference; `liftover` round-trips a BAM between
transcript and genome coordinates, splicing reads at exon boundaries. `-a`
takes the `prepare` annotation table; `-f` (a `.fai`) is required for `t2g`:

- `coralsnake liftover -d t2g` (default) — transcript BAM → genome BAM (splices
  reads at exon boundaries, inserts introns).
- `coralsnake liftover -d g2t` — genome BAM → transcript BAM (clips to exons,
  joins spliced reads contiguously on the transcript).
- `--table` mode converts a tab-separated sites table (instead of a BAM): `t2g`
  reads a gene column + 1-based transcript position and appends
  `GenomeChrom`/`GenomePos`; `g2t` reads chrom/position/strand and appends
  `Gene`/`GenePos`.

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
# or fast table mode, with the `prepare` annotation table (-o output)
coralsnake annotate -i sites.tsv -o out.tsv --annotation annotation.tsv
```

### `metagene` — exon-aware metagene profiling

Built on the high-performance `polars` + `ruranges` stack. Computes the
distribution of sites relative to gene regions (5'UTR, CDS, 3'UTR) and can emit
binned statistics and a publication-ready profile plot. The `-p` profile plot
needs the `plot` extra; the tabular outputs (`-o`, `-s`, `--export-profile`)
work without it. `--export-profile FILE` writes the machine-readable profile
matrix TSV (`feature_type`, `feature_midpoint`, `count_*`) for downstream tools.

```bash
# Using a built-in reference (GRCh38) or a custom GTF:
coralsnake metagene -i sites.tsv.gz -r GRCh38 -H -m 1,2,3 -w 5 \
                    -o output.tsv -s scores.tsv -p plot.png

coralsnake metagene -i sites.bed -g custom.gtf.gz -m 1,2,3 -w 5 \
                    -o output.tsv -s scores.tsv -p plot.png
```

Manage the built-in references (groups `human` / `mouse` are supported;
`reference list` shows sizes — the `metagene --list/--download/--export-*`
flags are deprecated aliases):

```bash
coralsnake reference list
coralsnake reference download GRCh38        # or: human / mouse / all
```

One reference download serves every tool — the tools accept the reference
**by name** and reuse the cached parquet (auto-downloading it if missing);
`reference export` converts it to text views for external tools:

```bash
coralsnake metagene -i sites.tsv -r GRCh38 ...              # uses the parquet directly
coralsnake annotate -i sites.tsv -g GRCh38 -f GRCh38 ...    # cached GTF + linked genome
coralsnake liftover -i tx.bam -o genome.bam -a GRCh38 -f GRCh38 --sort
coralsnake motif -i sites.tsv -f GRCh38 ...
coralsnake reference export GRCh38 --table ref_table.tsv --gtf ref.gtf  # explicit text views
```

Genome FASTAs are too large to ship in the release, so they are **linked**:
`reference genome <ref>` (or `reference download <ref> --with-genome`)
streams the verified upstream genome (Ensembl/UCSC), decompresses and
indexes it under `~/.cache/coralsnake/genomes/`, and cross-checks the FASTA
headers against the reference's contig names.

The reference parquets are served from this repo's fixed `data` release and
cached in `~/.cache/coralsnake/`; rebuild or update them with
[`scripts/build_references.py`](scripts/build_references.py)
(see [scripts/README.md](scripts/README.md)).

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

- `motif` — fetch the strand-aware genomic sequence around each site, padded
  with `N` (`-n left,right` for asymmetric padding).
- `coordinate` — rename chromosome names between reference coordinate systems
  (built-in UCSC↔Ensembl mappings for hg38/mm39, or a custom `-m` TSV).

Both migrated from the standalone `variant` package with unchanged naming and
output format.

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

Builds a DNA/RNA sequence logo from a set of motif sequences. The scoring
engine is pure numpy; only the figure renderer needs matplotlib (`plot`
extra). `--matrix FILE` exports the position × base score matrix as TSV
(rows = motif positions, columns = `A C G T U` + any other symbol present) —
pure numpy, so it works without the `plot` extra, and the figure is optional
when a matrix is requested.

```bash
coralsnake logo -m ACGT -m ACGG -m CCGT -o logo.png
# or with per-motif weights from a file (seq\tcount)
coralsnake logo -i motifs.tsv -o logo.svg
# or export just the position x base score matrix as TSV (no plotting needed)
coralsnake logo -i motifs.tsv --matrix logo_scores.tsv
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
