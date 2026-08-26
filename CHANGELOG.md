# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.215] - 2026-08-27

### Added
- **`gbam2tbam`** (`coralsnake/gbam2tbam.py`, new CLI command): inverse of
  `liftover`/`tbam2gbam` - remaps a genome-aligned BAM onto transcript
  references. Clips reads to exons, maps to 5'->3' transcript coordinates,
  joins spliced reads contiguously (adjacent exons), ref-skips skipped exonic
  gaps, and flips strand flags for `-` transcripts. Reads with internal
  insertions/deletions or M-blocks extending into introns are safely skipped
  rather than emitting an invalid alignment.
- **`tbam2gbam`** is now its own CLI command (alongside the `liftover` alias).
- **`annotate --annotation <table>`**: fast precomputed-table mode (subsumes
  `annot`). `annotate` is now the single annotation tool - GTF mode
  (region + gene/transcript/pos + variant effect) and table mode share one
  engine.
- **CLI help**: `annot`/`effect` moved to a `Deprecated (use annotate)` panel;
  `tbam2gbam`/`gbam2tbam` added to `Read Mapping`.

## [0.0.214] - 2026-08-27

### Added
- **`coralsnake annotate`** (`coralsnake/annotate.py`): the **single annotation
  tool** merging `annot` + `effect`. One fixed schema (`gene_id transcript_id
  transcript_pos region gene_pos transcript_strand mut_type transcript_motif
  coding_pos codon_ref aa_pos aa_ref distance2splice`) from one GTF-based
  engine: labels a bare site with gene/transcript/position + region (no FASTA
  needed) and - with a genome FASTA + ref/alt - the full variant effect; also
  distinguishes **intronic** from **intergenic**. (`--annotation` table mode
  was added in 0.0.215.)
- **Grouped, coloured CLI help**: `COMMAND_GROUPS` organises commands into
  named panels (`Read Mapping`, `Site & Variant Annotation`, `Genomic Analysis`,
  `Visualization`).

### Fixed
- **Bogus `Error: 0` panel on subcommand `--help`**: `click.exceptions.Exit`
  (subclass of `RuntimeError`, raised by `--help`/`--version`) was being caught
  by `CoralsnakeGroup.invoke` and rendered as a fake error. `Exit`/`Abort` now
  propagate untouched; only genuine `ValueError`/`RuntimeError` are converted.
  (`CoralsnakeGroup` also now subclasses `rich_click.RichGroup`, restoring the
  styled help.)
- **`motif` boundary off-by-one** (`motif.get_motif`): for sites near the end of
  a chromosome the motif was one base too short (right-overhang) and misplaced
  the site base (both-overhang). Rewritten with a single clamp-and-pad formula;
  verified against an independent property check across **864 cases** (all
  positions x lpad x rpad x strand) - the old code failed 604, the fix passes
  all. Interior / left-overhang behaviour is byte-identical to before.
- **`prepare --with-txpos` minus-strand span** (`Transcript.to_tsv`): emitted an
  inverted `start > end` span for `-`-strand multi-exon genes. `transcript_start`
  / `transcript_end` are now always the genomic bounding box
  `[min_exon_start+1, max_exon_end]` for both strands.
- **`Mlogo` NaN weights leaked NaN scores** (`logo.__init__`): the min-shift was
  recomputed from the original (NaN-containing) weights array, defeating the
  `nan_to_num`. The shift now uses the cleaned array.
- **`liftover --sort` crashed on absolute output paths**: the temporary sort
  file was derived from the output basename and broke when the path contained a
  directory. Now uses `tempfile.mkstemp` in the target directory (with cleanup).
- **`overlap.annotate_with_features`**: fixed undefined `input_strands` /
  `feature_strands` in the `by_strand` branch, and the `keep_all=False` column
  selection that always expected a `Name` column (now gated on `annot_name`).
- **Experimental tests** (`tests/experimental/*`) no longer crash at collection
  when the benchmark data is absent; they skip cleanly.
- Dropped unused imports (`annotate._to_rna`, `gtf.numpy`, `gtf.ensure_dir`).

## [0.0.213] - 2026-08-26

### Changed
- **Unified error handling**: library functions raise (`ValueError`/`RuntimeError`)
  instead of `sys.exit`; the CLI (via a custom `CoralsnakeGroup`) surfaces these
  as concise `click.ClickException` messages instead of tracebacks.
- **ruranges group-ID helper**: `utils.interval_groups()` now centralises the
  repeated string-label → `uint32` group mapping used by `annotation`,
  `map_to_local`, `overlap` and `gtf` (fully vectorized via
  `np.unique(return_inverse=True)`).
- **`effect` classifier**: `_refine_cds_effect` now detects **synonymous**
  (same-amino-acid) substitutions as `Silent`.
- **`download`**: replaced raw `print()`/`input()` + ANSI `Colors` with the
  rich `Console` (stderr) used elsewhere.
- Removed dead code (`print_analysis_summary` + its module-level logger
  side-effect, `gtf.py` demo `__main__`), de-duplicated the annotation-flattening
  block in `tbam2gbam` into `_flatten_annotation`, and fixed the `liftover` help
  text.
- Version → 0.0.213.

## [0.0.212] - 2026-08-26

### Added
- **`motif`, `coordinate`, `effect`** top-level commands: the standalone `variant`
  package fused into coralsnake, with the old buggy dependencies removed:
  - `motif` — strand-aware motif fetch (pyfaidx → **pysam.FastaFile**).
  - `coordinate` — chromosome-name mapping (urllib3 → **stdlib urllib**).
  - `effect` — variant effect classification as a pure-Python classifier built
    on coralsnake's GTF/CDS machinery (pyensembl+varcode removed).
- Flat modules `coralsnake/motif.py`, `coordinate.py`, `effect.py` (with
  `Site`/`Annot`, IUPAC/CODON_TABLE and effect ordering in `effect.py`).
- Tests: `test_variant` (10 tests) covering constants, motif, coordinate,
  effect.

### Changed
- Version → 0.0.212; the fused commands read/write gzip and `-` (stdin/stdout)
  via `xopen`; progress logs go to stderr so piped stdout stays clean.
- **Flattened the `metagene/` subpackage** into flat top-level modules
  (`gtf.py`, `io.py`, `annotation.py`, `map_to_local.py`, `overlap.py`,
  `plotting.py`, `download.py`, `config.py`); `metagene/utils.py` was merged
  into `utils.py`. `coralsnake.metagene.X` → `coralsnake.X`.
- Merged shared helpers: the three metagene CSV click converters → one
  `_parse_csv`; `logo._require_plotting` + `plotting._require_plotting` →
  `utils.require_plotting`; `gbam2tbam` annotation flattening de-duplicated.
- Removed dead code (`gtf.py` demo `__main__`, unused `print_analysis_summary`
  and its module-level logger side-effect); fixed the `liftover` help text.
- Docs: added `docs/variant.md`, updated architecture/DESIGN/README/CLI/API
  docs to reflect the fused top-level commands and the flattened layout.

## [0.0.211] - 2026-08-26

### Added
- **`coralsnake metagene`** subcommand: full migration of the standalone
  `metagene` package into coralsnake, exposing metagene profiling
  (load_sites / load_gtf / map_to_transcripts / normalize_positions /
  map_to_local / overlap stats) plus optional plotting.
- **`coralsnake logo`** subcommand: DNA/RNA sequence-logo plotting, migrated
  from the standalone `motiflogo` package.
- `coralsnake/metagene` subpackage with the complete migrated metagene API.
- Optional `plot` extra (`coralsnake[plot]`) for matplotlib-based
  visualizations; the base install stays light (lazy imports raise a helpful
  error if matplotlib is missing).
- Tests: `test_metagene`, `test_map_to_local`, `test_logo`,
  `test_metagene_perf` (20+ tests).

### Changed
- **Dependencies upgraded** (`uv.lock`): `bwamem 0.0.52 -> 0.0.56`,
  `ruranges 0.1.1 -> 0.2.7`; added `polars>=1.30.0` as a runtime dependency.
  Verified behavior is identical across the upgrade (full suite + byte-identical
  `annot`/mapping output).
- **Performance**: `map_to_transcripts` now uses a vectorized sort +
  `group_by().first()` instead of a `group_by().map_groups()` Python apply
  (~20× faster, identical output). `map_to_local` uses
  `ruranges.numpy.group_cumsum` for strand-aware cumulative offsets
  (~7× faster, identical output).
- `Mlogo` score-matrix builder rewritten with vectorized numpy (`bincount` +
  codepoint lookup); also fixes a `0·log2(0)` NaN edge case in the 2-bit path.
- Switched `ruranges.overlaps` calls to the modern `ruranges.numpy.overlaps`
  API.
- **Code cleanup**: removed all dead imports / unused variables across the
  package (`ruff` clean).

## [0.0.210] - 2026-08-26

### Fixed
- Liftover: avoid negative intron CIGAR when lifting over overlapping/adjacent
  exons; fallback for gene ID when exactly one transcript exists.
