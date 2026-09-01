---
layout: default
title: Variant Analysis
nav_order: 7
---

# Variant Analysis

Coralsnake provides two variant-analysis commands — `motif` and `coordinate` —
plus the unified `annotate` command for site/variant annotation. `motif` and
`coordinate` are a migration of the standalone
[`variant`](https://github.com/y9c/variant) toolkit with the same naming and
output format, but built on coralsnake's `pysam` + `ruranges` stack (the old
`pyfaidx` / `urllib3` dependencies are gone).

All of these commands read/write gzip and `-` (stdin/stdout) transparently.

## `motif`

Fetch a genomic motif centred on each site (strand-aware, padded with `N`).

```bash
coralsnake motif -i sites.tsv -o motifs.tsv -f genome.fa -n 2,3 -w -H
```

| Option | Description                                   |
| ------ | --------------------------------------------- |
| `-i, --input` | Input position file.                     |
| `-o, --output` | Output file.                            |
| `-f, --fasta` | Reference FASTA (**required**).           |
| `-n, --npad` | Padding bases; comma for left,right pads (`2,3`). |
| `-u, --to-upper` | Uppercase the motif.                   |
| `-w, --wrap-site` | Wrap the site in `[ ]`.                |
| `-H, --with-header` | Input has a header.                    |
| `-c, --columns` | Chrom,Pos,Strand columns (default `1,2,3`). |

## `coordinate`

Map chromosome names between reference coordinate systems (UCSC ↔ Ensembl,
custom, or built-in).

```bash
coralsnake coordinate -i sites.tsv -o mapped.tsv -m chrom_map.tsv -c 1
coralsnake coordinate -i sites.tsv -o mapped.tsv -M U2E
```

| Option | Description                                      |
| ------ | ------------------------------------------------ |
| `-m, --reference-mapping` | Custom 2-column mapping file.            |
| `-M, --buildin-mapping` | `U2E`, `E2U`, `U2E-hg38`, `E2U-hg38`, `U2E-mm39`, `E2U-mm39`. |
| `-c, --columns` | Chrom column (default `1`).                    |
| `-k, --keep-original` | Append the rename instead of replacing.   |
| `-H, --with-header` | Input has a header.                        |

## `annotate`

The unified site / variant annotation command — one GTF-based engine, one fixed
output schema. Use `annotate` for anything that falls under "which
gene/transcript is this, and what does it do?".

The same command serves a bare site (only chrom,pos,strand given ->
gene/transcript/position + region, no FASTA needed) **and** a full variant
(chrom,pos,strand,ref,alt + a genome FASTA -> coding codon/AA + effect):

```bash
# site labeling (region + gene/transcript/position) - just a GTF
coralsnake annotate -i sites.tsv -o out.tsv --reference-gtf annotation.gtf -c 1,2,3

# variant effect - add a genome FASTA (+ ref/alt columns)
coralsnake annotate -i variants.tsv -o effects.tsv \
                    --reference-gtf annotation.gtf --reference-transcript genome.fa -s -a
```

| Option | Description                              |
| ------ | ---------------------------------------- |
| `-i, --input` | Input site/variant file.          |
| `-o, --output` | Output annotation file.          |
| `-g, --reference-gtf` | Reference GTF (**required**).     |
| `-f, --reference-transcript` | Genome FASTA (needed for motif/codon/AA). |
| `-s, --strandness` | Use strand information.          |
| `-a, --all-effects` | Output all overlapping effects (not just the top). |
| `-n, --npad` | Padding bases (default 10).       |
| `-H, --with-header` | Input has a header.            |
| `-c, --columns` | Chrom,Pos,Strand,Ref,Alt (default `1,2,3,4,5`; Ref/Alt optional). |

### Output columns (single fixed schema)

```
gene_id  transcript_id  transcript_pos  region  gene_pos  transcript_strand
mut_type  transcript_motif  coding_pos  codon_ref  aa_pos  aa_ref  distance2splice
```

`region` (Intergenic / Intronic / 5'UTR / CDS / 3'UTR / Noncoding) is always
filled. The `mut_*`/coding columns are filled only when a genome FASTA and
ref/alt were supplied; otherwise they are empty.

## Python API

The commands are also importable functions.

```python
from coralsnake.annotate import run_annotate, Annotation  # unified
from coralsnake.motif import get_motif
from coralsnake.coordinate import run_coordinate
```
