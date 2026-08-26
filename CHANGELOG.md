# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.212] - 2026-08-26

### Added
- **`coralsnake variant`** command group: migration of the standalone `variant`
  package, with the old buggy dependencies removed:
  - `variant motif` — strand-aware motif fetch (pyfaidx → **pysam.FastaFile**).
  - `variant coordinate` — chromosome-name mapping (urllib3 → **stdlib urllib**).
  - `variant effect` — variant effect classification as a pure-Python
    classifier built on coralsnake's GTF/CDS machinery (pyensembl+varcode
    removed).
- `coralsnake/variant` subpackage (Site/Annot dataclasses, IUPAC/CODON_TABLE,
  effect ordering).
- Tests: `test_variant` (10 tests) covering constants, motif, coordinate,
  effect.

### Changed
- Version → 0.0.212; setuptools packages include `coralsnake.variant`.
- Docs: added `docs/variant.md`, updated architecture/DESIGN/README/CLI/API
  docs to reflect the migrated variant subcommands.

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
