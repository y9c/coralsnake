# Built-in reference data — build & update

`coralsnake metagene --download <ref>` fetches precomputed reference parquet
files (one per genome build) from the **`data` release of this repository**
(see `coralsnake/config.py`). They are cached under
`~/.cache/coralsnake/<REF>.parquet`.

Each parquet is the exon-level reference produced by
`coralsnake.gtf.prepare_exon_ref` (the same schema a custom GTF produces
locally), compressed with ZSTD.

## The `data` release is fixed

By design, the `data` release (currently a prerelease) is **immutable**: its
assets are never replaced. When the reference data needs an update (new
Ensembl/UCSC release, bugfix), publish under a **new tag** (e.g. `data-v2`)
and bump `GITHUB_DOWNLOAD_BASE` in `coralsnake/config.py` in the same change,
so released coralsnake versions always point at the data they were tested
against.

## Building the parquets

`build_references.py` is driven by `coralsnake.config.BUILTIN_REFERENCES`
(single source of truth: reference name, source GTF path, output filename).

```bash
# List the expected source files (and whether they are present)
python scripts/build_references.py --list

# Build everything (source GTFs under <source-dir>, layout = source_file;
# --fetch downloads any missing sources from the recorded source_url)
python scripts/build_references.py --all --fetch --source-dir ./sources --out data/parquet

# Build one
python scripts/build_references.py --single sources/Homo_sapiens/raw/Homo_sapiens.GRCh38.110.gtf.gz GRCh38
```

Every build runs a schema check (column set of `prepare_exon_ref`) and an
invariant check (exon offsets within transcript length, last exon ends at the
transcript end) before writing, so a silent logic regression cannot ship a
corrupt reference.

## Source GTFs

Place the raw GTFs under `<source-dir>` with the relative layout recorded in
`BUILTIN_REFERENCES[ref]["source_file"]` (e.g.
`Homo_sapiens/raw/Homo_sapiens.GRCh38.110.gtf.gz`) — or simply run the build
with `--fetch`: every reference carries a verified `source_url` in
`coralsnake/config.py`, and the script downloads whatever is missing.

Upstream locations (verified 2026-09-02; Ensembl is now at release 116):

- Ensembl (animals + *S. cerevisiae*):
  `https://ftp.ensembl.org/pub/release-<N>/gtf/<species>/<file>.gtf.gz`
- Ensembl Genomes (plants + *S. pombe*; `current` is frozen at release 63):
  `https://ftp.ebi.ac.uk/ensemblgenomes/pub/<plants|fungi>/current/gtf/<species>/<file>.gtf.gz`
- UCSC (gene GTFs moved under `bigZips/genes/`; the old top-level `genes/`
  dirs and dated combined GTFs are gone):
  `https://hgdownload.soe.ucsc.edu/goldenPath/<asm>/bigZips/genes/<file>.gtf.gz`

Note: v1's UCSC references were built from dated combined GTFs (e.g.
`hg38.20221028.gtf.gz`) that are no longer served; the v2 build uses the
per-track `knownGene` (GENCODE) GTF for vertebrates and the best available
curated set (refGene / ncbiRefSeq) for model organisms. Ensembl-based
references are rebuilt from the same releases as v1.

## Publishing an update

```bash
# 1. Build the new set (see above), into e.g. data/parquet-v2/
# 2. Bump GITHUB_DOWNLOAD_BASE in coralsnake/config.py (e.g. .../data-v2)
# 3. Publish under the new tag:
python scripts/build_references.py --all --source-dir ./sources --out data/parquet-v2 \
    --publish --tag data-v2 --title "Reference data v2"
```

`--publish` refuses to touch an existing release unless `--force` is given
(clobber assets on the existing tag — only use that to re-serve unchanged
content after an upload failure).
