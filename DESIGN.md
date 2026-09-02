# Coralsnake Library Design Plan

Coralsnake is a general NGS genomics toolkit grown from a single-purpose
two-color transcriptome mapper. This document is the **design plan** for
developers extending it. (A web-friendly version lives in
[`docs/architecture.md`](docs/architecture.md).)

## Goals

1. **Performance-first.** Hot paths must use vectorized primitives
   (ruranges Rust kernels, polars, numpy) — never per-interval Python loops.
2. **Thin CLI, reusable library.** Every subcommand is a thin wrapper around an
   importable function.
3. **Strand-aware, coordinate-correct.** Consistent 0-based half-open internal
   coordinates and a single strand convention.
4. **Light base install.** Heavy deps (matplotlib) are optional extras.
5. **Byte-identical output on refactor.** Behaviour-changing rewrites are
   verified against the prior output before landing.

## Module map

| Module | Responsibility | Key deps |
| ------ | -------------- | -------- |
| `cli.py` | command group, option parsing | rich-click |
| `seqops.c` | C kernels: RC, batch conversion, `score_and_tag`, `reverse_md` | (compiled) |
| `gtf2tx.py` | GTF/GFF → transcript FASTA | pysam |
| `tbam2gbam.py` / `gbam2tbam.py` | BAM liftover | pysam, ruranges |
| `genegroup.py` | gene clustering / consensus | numpy, scipy |
| `rnaseqc.py` | RNA-seq QC metrics (mapping + exon classification + counts) | pysam, ruranges, polars, xopen |
| `annot.py` | site → gene/transcript annotation | ruranges, polars-free I/O |
| `utils.py` | `Span`, `Transcript`, helpers | pysam, numpy |
| `logo.py` | sequence logo (numpy scoring) | numpy, (matplotlib opt.) |
| `motif.py` | genomic motif fetch (strand-aware) | pysam |
| `coordinate.py` | chrom-name mapping (UCSC↔Ensembl) | stdlib urllib |
| `effect.py` | variant-effect annotation + `Annot`/`Site`, genetic code | polars, pysam |
| `gtf.py` / `io.py` | metagene GTF parsing / site & reference I/O | polars, ruranges |
| `annotation.py` / `map_to_local.py` | metagene profiling core | polars, ruranges |
| `overlap.py` / `plotting.py` | metagene bin stats / plotting | polars, numpy |
| `download.py` / `config.py` | built-in reference download & config | stdlib urllib |

## Core data structures (`utils.py`)

- `Span(start, end)` — immutable half-open interval.
- `Transcript` — ordered exon collection with strand-aware helpers
  (`cum_exon_lens`, `length`, `get_seq`, `get_gene_spans`); the backbone for
  both `gtf2tx` and `tbam2gbam`.

## The ruranges contract

- Always use `from ruranges.numpy import ...` (modern API). The top-level
  `ruranges.overlaps` no longer exists in current releases.
- `group_cumsum(negative_strand=<bool>)` encodes **plus strand as True**.
- Kernels take 1-D numpy arrays of equal length; groups are `uint32` ids.

## Metagene pipeline

```
GTF + sites ─▶ map_to_transcripts ▶ normalize_positions ▶ gene_bins/stats/splits
```

- `map_to_transcripts`: overlap sites to exons, then pick the best transcript
  per gene via a *vectorized* sort + `group_by().first()` (NOT a python apply).
- `map_to_local`: strand-aware cumulative offsets via `ruranges.group_cumsum`.
- Caching: `load_gtf` writes a sidecar `.parquet`.

## Optional-dependency pattern

```python
def _require_plotting():
    try:
        import matplotlib
        return ...
    except ImportError as e:
        raise ImportError("Install with: pip install 'coralsnake[plot]'") from e
```

`coralsnake[plot]` is the only extra. All plotting is lazy.

## Testing strategy

- `tests/` mirrors `coralsnake/`; `tests/test_metagene_perf.py` guards the
  fast paths with generous time ceilings.
- Vectorized rewrites are validated by comparing output to the previous
  implementation (see `annotation.py` / `map_to_local.py` history).

## `variant` absorbed (fused into the package)

The standalone [`variant`](https://github.com/y9c/variant) toolkit has been
fused into coralsnake as flat top-level modules with its old deps removed:

| Module / command | Old dep (removed) | Replacement |
| ---------------- | ----------------- | ----------- |
| `motif.py` / `motif` | pyfaidx | pysam.FastaFile |
| `coordinate.py` / `coordinate` | urllib3 | stdlib urllib |
| `effect.py` / `effect` | pyensembl, varcode | pure Python + coralsnake GTF |

Naming and output format are unchanged; `effect.py` also carries the shared
`Annot`/`Site` dataclasses and the genetic-code / IUPAC tables. The classifier
reuses coralsnake's GTF/CDS machinery and, like the rest of the package,
streams I/O through `xopen` (gzip + `-` stdin/stdout).

Where coralsnake already provided a primitive (reverse_complement, GTF→exon/
codon indexing, gzip I/O), the migration fuses instead of duplicating.
