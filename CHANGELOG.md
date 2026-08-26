# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
