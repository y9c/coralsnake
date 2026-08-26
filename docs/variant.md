---
layout: default
title: Variant Analysis
nav_order: 7
---

# Variant Analysis

Coralsnake provides three variant-analysis commands — `motif`, `coordinate`
and `effect` — fused into the top-level CLI. They are a migration of the
standalone [`variant`](https://github.com/y9c/variant) toolkit with the same
naming and output format, but built on coralsnake's `pysam` + `ruranges` stack
(the old `pyfaidx` / `urllib3` / `pyensembl` + `varcode` dependencies are gone).

All three commands read/write gzip and `-` (stdin/stdout) transparently.

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

## `effect`

Annotate the effect of each variant (chrom, pos[, strand, ref, alt]). It
classifies the site into a region (5'UTR / CDS / 3'UTR / splice), computes
transcript-relative coordinates, and (for coding positions) the codon +
amino-acid context.

```bash
coralsnake effect -i variants.tsv -o effects.tsv \
                  --reference-gtf annotation.gtf \
                  --reference-transcript transcripts.fa \
                  -s -a
```
| Option | Description                              |
| ------ | ---------------------------------------- |
| `-i, --input` | Input variant file.                |
| `-o, --output` | Output annotation file.           |
| `-g, --reference-gtf` | Reference GTF (**required**).     |
| `--reference-transcript` | Transcript FASTA (repeatable). |
| `-s, --strandness` | Use strand information.          |
| `-a, --all-effects` | Output all effects (not just top). |
| `-u, --pU-mode` | Prioritise RNA genes.             |
| `-n, --npad` | Padding bases (default 10).         |
| `-H, --with-header` | Input has a header.            |
| `-c, --columns` | Chrom,Pos,Strand,Ref,Alt (default `1,2,3,4,5`). |

### Output columns

```
mut_type  gene_type  gene_name  gene_pos  transcript_name  transcript_pos
transcript_motif  transcript_strand  coding_pos  codon_ref  aa_pos  aa_ref
distance2splice
```

The field order is identical to the standalone `variant` package.

## Python API

The commands are also importable functions.

```python
from coralsnake.effect import Annot, Site, expand_base, reverse_base
from coralsnake.motif import get_motif
from coralsnake.coordinate import run_coordinate
from coralsnake.effect import build_transcript_index, run_effect
```
