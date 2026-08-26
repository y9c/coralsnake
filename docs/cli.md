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
| `prepare`  | Extract primary transcript from GTF/GFF file.          |
| `map`      | Map reads to reference genome using BWA-MEM.           |
| `liftover` | Remap transcriptome-aligned reads to genome BAM.       |
| `annot`    | Annotate a TSV of genomic sites with transcript info.  |
| `group`    | Group genes and build consensus sequences.             |
| `metagene` | Metagene profiling across 5'UTR/CDS/3'UTR.             |
| `logo`     | Plot a DNA/RNA sequence-logo.                          |
| `variant`  | Variant analysis: motif / coordinate / effect.         |

Run `coralsnake <command> --help` for the full option list of any command.

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

## `map`

Map paired/single-end reads to a reference with BWA-MEM, tuned for dual-base
conversion chemistry (MK/KM).

```
coralsnake map -1 R1.fq.gz -2 R2.fq.gz -r ref.fa -o out.bam \
               --fwd-ref --max-a2g-ratio 0.1 --threads 8
```

| Option | Description                                        |
| ------ | -------------------------------------------------- |
| `-1, --r1-file` | Read 1 FASTA/Q.                                 |
| `-2, --r2-file` | Read 2 FASTA/Q.                                 |
| `-r, --ref-file` | Reference FASTA (repeatable).                   |
| `-o, --output-file` | Output BAM (repeatable, matches refs).        |
| `-u, --unmap-file` | BAM for unmapped reads.                       |
| `--report` | HISAT2-style mapping summary (use `-` for stdout). |
| `-m, --max-mismatches` | Max bad mismatches (default 10).            |
| `-t, --threads` | Worker processes (default 8).                 |
| `--min-alignment-length` | Default 20.                             |
| `--min-mapping-ratio` | Default 0.5.                              |
| `--max-a2g-ratio` | Max A→G proportion (default 1.0).             |
| `--max-c2t-ratio` | Max C→T proportion (default 1.0).             |
| `--index-dir` | BWA index directory (repeatable).              |
| `--index-only` | Build indices without mapping.                |
| `--fwd-lib/--rev-lib` | Library orientation.                     |
| `--fwd-ref/--rev-ref/--dbl-ref` | Reference strand.                  |

---

## `liftover`

Remap a transcriptome-aligned BAM back to genome coordinates.

```
coralsnake liftover -i tx.bam -o genome.bam -a annot.tsv -f faidx.fai --sort
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-i, --input-bam` | Input BAM (**required**).               |
| `-o, --output-bam` | Output BAM (**required**).              |
| `-a, --annotation-file` | Annotation TSV (**required**).         |
| `-f, --faidx-file` | Reference `.fai` (**required**).          |
| `-t, --threads` | Threads.                                   |
| `-s, --sort` | Sort output by coordinate.                    |

---

## `annot`

Annotate a TSV of genomic sites (chrom, position, strand) with the matching
gene/transcript and the transcript-relative position.

```
coralsnake annot -i sites.tsv -o annotated.tsv -a annot_table.tsv \
                 -c 1,2,3 -k
```

| Option | Description                                  |
| ------ | -------------------------------------------- |
| `-i, --input-file` | Input TSV (**required**).               |
| `-o, --output-file` | Output TSV (**required**).              |
| `-a, --annot-file` | Annotation table (**required**).          |
| `-c, --cols` | Columns of Chrom,Pos,Strand (default `1,2,3`). |
| `-k, --keep-na` | Emit `.` for unannotated sites.           |
| `-l, --collapse-annot` | Collapse multiple hits into one row.   |
| `-n, --add-count` | Append a hit count column.               |
| `-H, --skip-header` | Input has a header line.                 |

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

## `variant`

Genomic variant analysis, migrated from the standalone `variant` package.
Three subcommands — `motif`, `coordinate` and `effect` — keep the original
naming and output format but are built on coralsnake's `pysam` + `ruranges`
stack (the old `pyfaidx` / `urllib3` / `pyensembl` + `varcode` dependencies
are gone, see [Design](architecture.md)).

```
usage: coralsnake variant COMMAND [OPTIONS]
```

### `variant motif`

Fetch a genomic motif (reference sequence) around each input site, with
strand-aware reverse-complementing and `N`-padding at chromosome ends.

```
coralsnake variant motif -i sites.tsv -o motifs.tsv -f genome.fa \
                         -n 10 -c 1,2,3 -u -w
```

| Option | Description                                      |
| ------ | ------------------------------------------------ |
| `-i, --input` | Input position file (default `-`).   |
| `-o, --output` | Output file (default `-`).          |
| `-f, --fasta` | Reference FASTA (**required**).          |
| `-n, --npad` | Padding bases (default `10`; use `2,3` for L,R). |
| `-H, --with-header` | Input has a header line.            |
| `-c, --columns` | Site columns Chrom,Pos,Strand (default `1,2,3`). |
| `-u, --to-upper` | Emit motif in upper case.            |
| `-w, --wrap-site` | Wrap the motif site in `[...]`.      |

### `variant coordinate`

Map chromosome names between coordinate systems (UCSC ↔ Ensembl). Accepts a
custom 2-column mapping file or a built-in preset.

```
coralsnake variant coordinate -i in.tsv -o out.tsv -M U2E -c 1 -k
```

| Option | Description                                      |
| ------ | ------------------------------------------------ |
| `-i, --input` | Input file (default `-`).              |
| `-o, --output` | Output file (default `-`).           |
| `-m, --reference-mapping` | 2-col mapping file `chrom<TAB>renamed`. |
| `-M, --buildin-mapping` | `U2E`, `E2U`, `U2E-hg38`, `E2U-hg38`, `U2E-mm39`, `E2U-mm39`. |
| `-c, --columns` | Chrom column (default `1`).            |
| `-H, --with-header` | Input has a header line.            |
| `-k, --keep-original` | Keep original chrom as an extra column. |

### `variant effect`

Annotate each variant with its predicted effect (5'UTR/CDS/3'UTR, splice
sites, intronic, intergenic, …), the affected gene/transcript, codon and
amino-acid context. Uses a pure-Python classifier over coralsnake's GTF
machinery — no `pyensembl`/`varcode`.

```
coralsnake variant effect -i sites.tsv -o effects.tsv \
                          -g annotation.gtf -f genome.fa -s -n 10 -a
```

| Option | Description                                      |
| ------ | ------------------------------------------------ |
| `-i, --input` | Input file (default `-`).              |
| `-o, --output` | Output annotation file (default `-`). |
| `--reference-gtf` | Reference GTF (**required**).        |
| `--reference-transcript` | Reference transcript FASTA(s).  |
| `--reference-protein` | Reference protein FASTA(s).       |
| `-e, --release` | Ensembl release (for built-ins).     |
| `-s, --strandness` | Use strand information.             |
| `-u, --pU-mode` | Prioritise RNA genes.                |
| `-n, --npad` | Padding bases for motif (default 10).          |
| `-a, --all-effects` | Output all effects (not just the top). |
| `-H, --with-header` | Input has a header line.            |
| `-c, --columns` | Site columns Chrom,Pos,Strand,Ref,Alt (default `1,2,3,4,5`). |
