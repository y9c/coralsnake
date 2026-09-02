---
layout: default
title: CLI Reference
nav_order: 3
---

# CLI Reference

```
Usage: coralsnake [OPTIONS] COMMAND [ARGS]...
```

## Commands

| Command    | Description                                            |
| ---------- | ------------------------------------------------------ |
| `reference`| Manage built-in references (list/download/export/genome). |
| `prepare`  | Extract primary transcript from GTF/GFF file.          |
| `liftover` | Bidirectional BAM liftover (`--direction t2g` default / `g2t`). |
| `annotate` | Unified site/variant annotation (region + gene/effect).|
| `group`    | Group genes and build consensus sequences.             |
| `metagene` | Metagene profiling across 5'UTR/CDS/3'UTR.             |
| `logo`     | Plot a DNA/RNA sequence-logo.                          |
| `motif`    | Fetch a genomic motif around variant sites.            |
| `coordinate` | Map chromosome names between coordinate systems.    |

Run `coralsnake <command> --help` for the full option list of any command.

> **Mapping is out of scope:** `coralsnake map` was removed. Align reads with
> `bwamem map` (plain) or `prismalign map` (nucleotide-conversion / two-color),
> then `prepare` + `liftover` for the exon-aware steps.

---

## `reference`

Manage the built-in exon references (exon-level parquets built from canonical
GTFs). **One download serves every tool**: `metagene -r`, `liftover
-a/--annotation-file`, `annotate --annotation` / `--reference-gtf`,
`motif -f` and `annotate -f` all accept a built-in reference **name** (e.g.
`GRCh38`) and reuse the cached object, auto-downloading the parquet and
deriving/caching the table or GTF on first use.

```
coralsnake reference list
coralsnake reference download GRCh38            # or: human / mouse / all
coralsnake reference download GRCh38 --with-genome
coralsnake reference genome GRCh38
coralsnake reference export GRCh38 --table ref_table.tsv --gtf ref.gtf
```

| Subcommand | Description |
| ---------- | ----------- |
| `list` | List the built-in references with their download sizes. |
| `download <ref\|human\|mouse\|all>` | Download reference parquet(s) into the cache (atomic writes, streaming progress). `--with-genome` also fetches the linked genome FASTA(s). |
| `genome <ref>` | Download the linked genome FASTA (+ `.fa.fai` index) into `~/.cache/coralsnake/genomes/<ref>.fa`. Genome sequences are too large to ship in the `data` release (human ≈ 900 MB compressed), so each reference records a verified upstream URL (`config.GENOME_URLS`); the FASTA is decompressed, indexed and cross-checked against the reference's contig names (hard error on zero overlap). |
| `export <ref> --table FILE --gtf FILE` | Export the cached reference as the `prepare` annotation table (feeds `liftover -a` / `liftover --table` / `annotate --annotation`) and/or a GTF (feeds `annotate --reference-gtf`), for use by external tools. |

---

## `prepare`

Extract the primary transcript from a GTF/GFF file into a reference
transcript FASTA.

```
coralsnake prepare -g annotation.gtf -f genome.fa -o transcripts.fa \
                   --with-codon --with-genename --filter-biotype protein_coding
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-g, --gtf-file` | GTF file (**required**).                |
| `-f, --fasta-file` | Reference FASTA.                        |
| `-o, --output-file` | Output FASTA (**required**).           |
| `-s, --seq-file` | Sequence file (requires `--fasta-file`).  |
| `-c, --with-codon` | Include codon info.                    |
| `-n, --with-genename` | Include gene name.                   |
| `-t, --with-biotype` | Include biotype.                     |
| `-x, --with-txpos` | Include transcript position.            |
| `-z, --sanitize` | Sanitize sequence names.                |
| `-U/-u, --seq-upper/--seq-lower` | Case of sequence.       |
| `-l, --line-length` | Output line length (0 = no wrap).      |

---

## `refine`

Clean a genome FASTA and/or GTF **before** `prepare`: seqname renaming/filtering,
gene/transcript name and type normalization, missing gene/transcript/exon row
creation, overlapping-exon merge, canonical-transcript flagging, and coordinate
checks. Codon and UTR features are preserved (so the output stays usable for
`metagene` and `annotate`). Indexing uses pysam — no external samtools/bgzip/tabix.

The GTF side is built on the shared `coralsnake.genemodel.GeneModel`
object — the same read/serialize layer `prepare` uses — so a refined GTF is a
drop-in replacement for the input:

- **Biotypes written for `prepare`**: both `gene_type`/`transcript_type`
  (GENCODE style) and `gene_biotype`/`transcript_biotype` (Ensembl style) are
  emitted, so `prepare --filter-biotype / --with-biotype` work on the refined
  output even when the input only carried `*_type`.
- **Canonical transcripts are flagged** with `is_canonical` *and* tagged
  `Ensembl_canonical`, so `prepare`'s ranking (MANE/Ensembl-canonical priority)
  picks the transcript `refine` selected (longest protein-coding, or from
  `--canonical-transcripts`).

```
coralsnake refine -f genome.fa -g annotation.gtf -o outdir -n hg38 \
                  -m chrom_map.tsv -c canonical.tsv
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-f, --fasta-file` | Input genome FASTA (`.gz` ok).        |
| `-g, --gtf-file` | Input GTF (`.gz` ok).                     |
| `-o, --outdir` | Output directory (default `./`).           |
| `-n, --name` | Output name prefix (default: outdir name). |
| `-m, --rename-mapper` | TSV: old seqname → new seqname.    |
| `-p, --seqname-pattern` | Keep seqnames matching this regex.    |
| `-c, --canonical-transcripts` | TSV of canonical transcript IDs (1st column). |

Outputs (under `<outdir>/<name>.*`): `.genome.fasta` (+ `.fai`, `.genome.sizes`),
`.annotation.gtf` (+ `.gz` + tabix index), `.skip.gtf` (genes that failed
checks), `.gene_features_summary.txt`. At least one of `-f`/`-g` is required.

---


## `liftover`

Remap a BAM between transcript and genome coordinates, **both directions**:

- `--direction t2g` (default) — transcriptome-aligned BAM → genome coordinates.
- `--direction g2t` — genome-aligned BAM → transcript coordinates (the former
  standalone `gbam2tbam` command, now fully fused into `liftover`).

```
coralsnake liftover -i tx.bam -o genome.bam -a annot.tsv -f faidx.fai --sort
coralsnake liftover -d g2t -i genome.bam -o tx.bam -a annot.tsv --sort
coralsnake liftover -i tx.bam -o genome.bam -a GRCh38 -f GRCh38 --sort
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-d, --direction` | `t2g` (default) or `g2t`.             |
| `-i, --input-bam` | Input BAM (**required**).               |
| `-o, --output-bam` | Output BAM (**required**).              |
| `-a, --annotation-file` | Annotation TSV, or a built-in reference name (cached table) (**required**). |
| `-f, --faidx-file` | Reference `.fai`, or a built-in reference name (required for `t2g`). |
| `-t, --threads` | Threads.                                   |
| `-s, --sort` | Sort output by coordinate.                    |

---

## `annotate`

The unified site / variant annotation command, with a single output schema.
One command serves a bare site (chrom,pos,strand → gene/transcript/position +
region) and a full variant (chrom,pos,strand,ref,alt + genome FASTA → codon/AA +
effect). A precomputed `prepare`-style table (`--annotation`) gives fast table
mode.

```
coralsnake annotate -i sites.tsv -o out.tsv --reference-gtf annotation.gtf -c 1,2,3
coralsnake annotate -i variants.tsv -o effects.tsv \
                    --reference-gtf annotation.gtf \
                    --reference-transcript genome.fa -s -a
coralsnake annotate -i sites.tsv -o out.tsv --annotation prepared_table.tsv
coralsnake annotate -i variants.tsv -o effects.tsv -g GRCh38 -f GRCh38 -s -a
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-i, --input` | Input site/variant file (default `-` = stdin; **1-based** positions, like `motif`). |
| `-o, --output` | Output file (default `-` = stdout).       |
| `-g, --reference-gtf` | Reference GTF, or a built-in reference name (cached GTF) (GTF mode). |
| `--annotation` | Precomputed `prepare` table, or a built-in reference name (fast table mode). |
| `-f, --reference-transcript` | Genome FASTA, or a built-in reference name (linked genome, fetched on demand) — needed for motif/codon/AA. |
| `-s, --strandness` | Use strand information.                |
| `-n, --npad` | Padding bases for motif (default 10).        |
| `-a, --all-effects` | Output all overlapping effects (not just top). |
| `-H, --with-header` | Input has a header line.            |
| `-c, --columns` | Chrom,Pos,Strand,Ref,Alt (default `1,2,3,4,5`; Ref/Alt optional). |

Output schema: `gene_id transcript_id transcript_pos region gene_pos
transcript_strand mut_type transcript_motif coding_pos codon_ref aa_pos aa_ref
distance2splice`. The `mut_*`/coding columns are filled only when a genome
FASTA and ref/alt were supplied.
---

## `metagene`

Run metagene profiling: the distribution of genomic sites relative to gene
regions (5'UTR / CDS / 3'UTR), with optional binned score statistics and a
publication-ready profile plot (requires `coralsnake[plot]`). Use a built-in
reference (`-r GRCh38`) or a custom GTF (`-g file.gtf`). Reference
management lives under the top-level [`reference`](#reference) command; the
`--list` / `--download` / `--export-table` / `--export-gtf` flags on
`metagene` are deprecated aliases kept for backward compatibility.

```
coralsnake metagene -i sites.tsv.gz -r GRCh38 -H -m 1,2,3 -w 5 \
                    -o out.tsv -s scores.tsv -p plot.png
coralsnake metagene -i sites.bed -g custom.gtf.gz -m 1,2,3,6 -o out.tsv
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-i, --input` | Input file (BED/TSV/CSV).             |
| `-o, --output` | Output score table path.             |
| `-s, --output-score` | Output binned score statistics path. |
| `-p, --output-figure` | Metagene plot path (`coralsnake[plot]`). |
| `-r, --reference` | Built-in reference (e.g. `GRCh38`, `GRCm39`). |
| `-g, --gtf` | Custom GTF/GFF reference.          |
| `--region` | `all` (default) / `5utr` / `cds` / `3utr`. |
| `-b, --bins` | Number of bins (default 100).      |
| `-H, --with-header` | Input has a header line.     |
| `-S, --separator` | Input separator (default tab). |
| `-m, --meta-columns` | Coordinate columns (default `1,2,3,6`). |
| `-w, --weight-columns` | Weight/score columns (1-based). |
| `-n, --weight-names` | Names for the weight columns.   |
| `--score-transform` | `none` (default) / `log2` / `log10`. |
| `--normalize` | Normalize scores by transcript length. |
| `--list` | Deprecated: use `reference list`. List built-in references and exit. |
| `--download` | Deprecated: use `reference download`. Download a reference (`<name>`), a group (`human`, `mouse`), or `all`. |
| `--export-table FILE` | Deprecated: use `reference export --table`. Export the reference (with `-r`) as the `prepare` annotation table (TSV) and exit. |
| `--export-gtf FILE` | Deprecated: use `reference export --gtf`. Export the reference (with `-r`) as a GTF and exit. |

---

## `logo`

Plot a DNA/RNA sequence-logo from a set of motif sequences. The scoring engine
is pure numpy; rendering requires matplotlib (`pip install "coralsnake[plot]"`).

```
coralsnake logo -m ACGT -m ACGG -m CCGT -o logo.png
coralsnake logo -i motifs.tsv -o logo.svg   # one motif per line; seq<TAB>count for weights
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-m, --motifs` | Motif sequences (repeatable / comma-separated). |
| `-i, --input` | File of motifs, one per line.         |
| `-o, --output` | Output image (e.g. `logo.png`, `logo.svg`) (**required**). |
| `-w, --weights` | Comma-separated per-motif weights.  |
| `--t2u/--no-t2u` | Convert T→U (default on).        |
| `--2bit/--no-2bit` | 2-bit logo (default on).       |
| `--normed` | Normalize letter heights.               |

---

## `group`

Group related genes and build consensus sequences.

```
coralsnake group -f genes.fa -g genes.gtf -o grouped.tsv \
                 --output-consensus consensus.fa --threads 8
```

| Option | Description                                      |
| ------ | ------------------------------------------------ |
| `-f, --fasta-file` | FASTA file (**required**).                     |
| `-g, --gtf-file` | GTF file (**required**).                          |
| `-o, --output-file` | Grouped output TSV.                           |
| `-c, --output-consensus` | Consensus FASTA output.                   |
| `-r, --gene-name-regex` | Gene name regex.                           |
| `-b, --gene-biotype-list` | Biotype filter.                         |
| `-l, --gene-length-limit` | Max gene length (default 300).            |
| `-s, --cluster-threshold` | Clustering threshold (default 0.1).        |
| `-t, --threads` | Threads.                                        |

---

## `motif`

Fetch a genomic motif (reference sequence) around each input site, strand-aware
with `N`-padding at chromosome ends. Migrated from the standalone `variant`
package; the old `pyfaidx` dependency is gone (uses `pysam.FastaFile`).

```
coralsnake motif -i sites.tsv -o motifs.tsv -f genome.fa -n 10 -c 1,2,3 -u -w
coralsnake motif -i sites.tsv -o motifs.tsv -f GRCh38 -n 10 -c 1,2,3
```

| Option | Description                                      |
| ------ | ------------------------------------------------ |
| `-i, --input` | Input position file (default `-` = stdin).|
| `-o, --output` | Output file (default `-` = stdout).     |
| `-f, --fasta` | Reference FASTA, or a built-in reference name (linked genome, fetched on demand) (**required**). |
| `-n, --npad` | Padding bases (default `10`; use `2,3` for L,R). |
| `-H, --with-header` | Input has a header line.            |
| `-c, --columns` | Site columns Chrom,Pos,Strand (default `1,2,3`). |
| `-u, --to-upper` | Emit motif in upper case.            |
| `-w, --wrap-site` | Wrap the motif site in `[...]`.      |

---

## `coordinate`

Map chromosome names between coordinate systems (UCSC ↔ Ensembl). Accepts a
custom 2-column mapping file or a built-in preset.

```
coralsnake coordinate -i in.tsv -o out.tsv -M U2E -c 1 -k
```

| Option | Description                                      |
| ------ | ------------------------------------------------ |
| `-i, --input` | Input file (default `-` = stdin).       |
| `-o, --output` | Output file (default `-` = stdout).    |
| `-m, --reference-mapping` | 2-col mapping file `chrom<TAB>renamed`. |
| `-M, --buildin-mapping` | `U2E`, `E2U`, `U2E-hg38`, `E2U-hg38`, `U2E-mm39`, `E2U-mm39`. |
| `-c, --columns` | Chrom column (default `1`).            |
| `-H, --with-header` | Input has a header line.            |
| `-k, --keep-original` | Keep original chrom as an extra column. |
