# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **`scripts/build_references.py`** — builds the built-in metagene reference
  parquets from source GTFs (driven by `coralsnake.config.BUILTIN_REFERENCES`,
  using `prepare_exon_ref` + ZSTD) and can publish them to a GitHub data
  release (`--publish`). See `scripts/README.md` for the fixed-release
  convention and source GTF locations.

### Changed
- **Built-in reference data now hosted in this repo**: `coralsnake metagene
  --download` fetches from the `data` release of `y9c/coralsnake` instead of
  the no-longer-maintained `y9c/metagene` repo. File contents are byte-identical
  to the previous release (SHA-256 verified at migration).

### Added
- **`logo --matrix FILE`**: export the position × base score matrix as a TSV
  (rows = motif positions, columns = `A C G T U` + any other symbol present,
  values = per-base score from `Mlogo.scores`). Machine-readable interchange for
  downstream tools/report renderers; the figure is now optional, and matplotlib
  is no longer required when only the matrix is requested.
- **`metagene --export-profile FILE`**: write the metagene profile matrix TSV
  (`feature_type`, `feature_midpoint`, `count_*`) for downstream tools; the
  profile is computed whenever any of `-o`, `-s`, `-p` or `--export-profile` is
  requested (previously `-o` was mandatory).

## [Unreleased] (cont.) — `refine` command (from PR #2, reworked)

New `coralsnake refine` command (genome FASTA + GTF cleaning before `prepare`),
contributed by Zonggui Chen and reworked onto the package infrastructure:

- **Attribute parsing reuses the hardened `gtf2tx.parse_gtf_annot`** (handles
  quoted `;` inside values and a missing trailing `;`) instead of a hand-rolled
  parser.- **Codon/UTR features are preserved** (the PR dropped start/stop_codon and UTR
  rows, which would have broken `metagene` gene splits and `annotate` CDS
  effects).
- **Indexing via pysam** (`faidx`, `tabix_compress`/`tabix_index`) — no external
  samtools/bgzip/tabix on PATH, no `shell=True`.
- Crashes fixed: exon-only genes, CDS without a containing exon, and the empty
  default `name` (which produced hidden dotfile outputs).
- `refine` validates that at least one of `--fasta-file`/`--gtf-file` is given.
- Verified end-to-end: refining the yeast R64 GTF keeps all 145/145 codons and
  feeds `prepare --with-codon` correctly.

### New: shared `GeneModel` object (`coralsnake.genemodel`)

One object models one annotated genome (GTF/GFF3-style attributes): every
feature row is parsed once, kept losslessly, and grouped under its gene and
transcript (1-based GTF coordinates + 0-based span helpers). It is the single
source of truth for how coralsnake reads and writes gene annotations:

- `gtf2tx.read_gtf` (`prepare`) now consumes `GeneModel.iter_rows()` — the
  attribute regexes, comment handling and malformed-line rules live in one
  place; per-row ranking / biotype / GFF-Parent decisions are unchanged
  (`genegroup` and `prepare` output are byte-identical).
- `refine` is rebuilt on the object (load → mutate → `write_gtf` with sort +
  bgzip/tabix). The metagene/annotate path keeps its fast cached Polars loader.

Compatibility guarantees for the built-in tools:

- **Biotype aliases**: refined GTFs emit both `gene_type`/`transcript_type`
  (GENCODE) and `gene_biotype`/`transcript_biotype` (Ensembl), so
  `prepare --with-biotype / --filter-biotype` work on GENCODE-style inputs.
- **Canonical transcripts are tagged** `Ensembl_canonical` (plus `is_canonical`),
  so `prepare`'s ranking picks the transcript `refine` selected.
- Verified numerically identical for metagene/annotate: `prepare_exon_ref` on the
  original R64 GTF vs the refined one agrees on all 218 exon rows (coordinates,
  exon numbers, transcript offsets, codon positions, levels).

## [Unreleased]

### ⚠ Behavior change: `annotate` / `effect` input positions are now 1-based
- The fused `annotate` (GTF mode) and legacy `effect` commands interpreted the
  input `Pos` as **0-based**, while `motif`, `annot`, metagene 3-column mode,
  and `annotate --annotation` (table mode) all use the package-wide **1-based**
  site convention. `annotate`/`effect` now convert 1-based input positions at
  the boundary, so all site inputs behave identically. (This also fixes every
  coordinate-dependent output of the GTF mode — region, codon, motif, distance
  to splice — being shifted one base downstream for real inputs.)

### Performance
- **`annotate` (GTF mode) ~95× faster**: the per-site overlap used a Python
  bisect + linear scan over every gene-body span on the chromosome; it now
  uses chunked (50k-site) `ruranges` batch overlap, with the classifier running
  only on true overlaps. Output schema/order/tie-break unchanged.
- **`motif` batches FASTA access**: sites are grouped per chromosome and read
  through a sorted sliding 4 Mb window instead of one faidx seek per site
  (large win on big genomes); output order preserved.

### Bug fixes (all regression-tested)
- **`liftover -d g2t`**: reads on `-` transcripts now have SEQ/QUAL
  reverse-complemented (SEQ is stored reference-forward per the SAM spec — the
  previous output mismatched the transcript reference base-by-base);
  cross-chromosome transcript assignment fixed (a read on chrX could be mapped
  onto a chr1 transcript whose coordinates overlapped); mate position is
  remapped for paired reads when possible; CIGAR `P` op handled; stale `SA`/`XC`
  tags dropped.
- **`liftover -d t2g`**: CIGAR `=`/`X` (and `P`) ops now consume the reference
  (intron N's were missing → coordinates corrupted); `transcript_pos == length`
  boundary raises the intended `ValueError` and malformed reads are demoted to
  unmapped instead of aborting the run; mate position remapped (exact for `+`,
  approximated within one read length for `-`); `SA`/`XA` tags dropped.
- **`prepare`**: `--with-codon` now works for standard (Ensembl/GENCODE) GTFs
  whose codon lines lack `exon_number`; exon spans in the output table are now
  always in 5'→3' order (previously they flipped with `--seq-file`, corrupting
  minus-strand `transcript_pos` in `annot`/table mode).
- **metagene/GTF engine**: exon local offsets are correct even when a GTF lists
  exons out of transcript order (was: wrong/negative `Start_exon` for `-`
  genes); GTFs without codon features or `exon_number` no longer crash;
  `--weight-names` now tracks the renamed weight columns (was: profile computed
  on the wrong column); default `-m 1,2,3,6` auto-falls back to `1,2,3` for
  3-column site files; friendly errors for out-of-range `-m`/`-w` and `--bins
  <= 0`; warning when no sites have CDS coordinates.
- **`annotate`**: the stop codon's 3 bases are CDS (were partially 3'UTR, with
  the codon itself dropped from the CDS slice); placeholder `.`/`N` ref/alt no
  longer fabricate variant effects; equal-length multi-base changes are
  `ComplexSubstitution` (not `InFrameIndel`); untranslatable codons (gaps) are
  no longer mislabeled `Silent`; table mode emits uniform 16-field rows, no
  longer leaks the first data row into the header, and skips malformed rows;
  `--columns` lacking chrom/pos errors up front; stale `.pickle` annotation
  caches are re-parsed (mtime check); `effect` handles empty input / malformed
  rows; the deprecated `annot`/`effect` commands now print a runtime warning.
- **`io.load_sites`**: inputs with a header column named exactly
  `Chromosome`/`Start`/`End`/`Strand` no longer crash; `map_to_local` empty
  results keep the full column schema; reference downloads are atomic
  (temp file + rename, no more corrupted caches).
- **`motif`**: a site beyond the contig end yields an all-`N` motif instead of
  crashing the run; unknown chromosome names raise a clear error.
- **`overlap`**: duplicate input rows (same coordinates, different weights) are
  no longer collapsed; `type_ratios` length is validated.
- **`logo`**: non-ASCII motifs raise a clear error; import order fixed (ruff
  E402).
- **`seqops.c`**: `batch_base_conversion` validates list items (non-str items
  raised `SystemError`).
- **Verified NOT bugs** (audited, locked in by tests): minus-strand
  sequence/CIGAR/MD handling in `t2g` (matches the SAM convention, confirmed by
  samtools mpileup), the g2t reference walk direction, the gap-gap mask in
  `genegroup` (pyfamsa yields int arrays), the weighted-logo IC formula
  (inherited from `motiflogo`).

## [0.0.222] - 2026-08-31

### Architecture: `map` removed — N-color mapping lives in `prismalign`
- **`coralsnake map` is removed.** The two-color (and N-color) nucleotide-
  conversion mapping engine has moved out of both coralsnake *and* bwamem into
  a new dedicated package, **`prismalign`** (2-color MK/KM, 3-nt BS/SLAM/A2G and
  custom N-channel schemes), which uses **pluggable alignment backends** —
  `bwamem` (BWA-MEM), a pure-Python reference backend, and `mappy` (minimap2).
  bwamem itself is kept **light** (a thin BWA-MEM binding; its earlier color
  addition was reverted).
- **identity**: coralsnake is now the exon-aware RNA pipeline — `prepare`,
  `liftover` (both directions; splicing/joining), `annotate`, `metagene`,
  `variant`. README / DESIGN / metadata updated.

### CLI: `gbam2tbam` fused into `liftover`
- **`coralsnake gbam2tbam` is removed** — the g2t direction is now
  `coralsnake liftover --direction g2t`. `gbam2tbam.py` stays as the internal
  g2t implementation; `liftover` is the single BAM-conversion command
  (`-d t2g` transcript→genome, `-d g2t` genome→transcript).
- **Docs cleaned up**: the deprecated `annot`/`effect` commands are no longer
  documented (`annotate` is the one annotation tool); CLI/API/architecture
  pages now match the current command set, and an ASCII overview diagram was
  added to README + docs index.

### Speed (C)
- **`seqops.c` gains `reverse_md`** — a C kernel reversing MD:Z tags on `-`
  strand reads during `tbam2gbam` (pure-Python regex fallback kept). Tested 1:1
  against the Python implementation across deletion/mismatch/merge cases.
- **removed dead `mk_conversion`/`km_conversion` helpers** from `utils.py`
  (they called a non-existent `seqops.base_conversion`; unused).

## [0.0.221] - 2026-08-28

### Speed
- **`annotate` variant-effect path now caches each transcript's spliced sequence
  per run** (previously every site re-fetched/re-assembled the FASTA sequence of
  its transcript). With many sites in one gene this is a large speedup
  (benchmark: 20,000 sites in a single gene + FASTA in ~0.16 s).

### Stability
- **`annotate` skips malformed input rows instead of aborting**: a row with a
  missing chromosome/position column or a non-integer position is now skipped
  and processing continues (previously it aborted the whole run with a
  `ValueError`).

## [0.0.220] - 2026-08-28

### Added / fixed
- **`gbam2tbam` now handles reference deletions (`D`)** (the last previously-
  skipped case): exonic `D` bases are mapped to transcript positions and emitted
  as `D` on the transcript. Added an exact `-`-strand orientation proof (a real
  reverse-on-genome read, flag 0x10, comes out forward on the transcript, 0x00,
  at the correct coords). Also added exact internal-soft-clip/insertion and
  intron-dipping tests.
- **Randomized property/fuzz tests** (`tests/test_properties.py`): seeded
  invariants over thousands of cases for `motif` (length/center/content vs an
  independent reference), `gbam2tbam.remap_read` (any produced read is
  structurally valid: query length == M+I+S, in-bounds) on both strands, and
  `_assemble_cigar` (length + reference-span invariants).
- **Verified the built-in-reference download URLs resolve** (they point at the
  `y9c/metagene` release assets, which serve the parquet files correctly - no
  change needed).

## [0.0.219] - 2026-08-27

### Fixed / improved
- **`gbam2tbam` read-completeness**: previously it silently dropped any read
  whose CIGAR had a soft-clip/insertion in the middle, or whose M blocks ran
  into an intron. It now walks the full CIGAR and converts these correctly:
  internal soft-clips become `S`, insertions become `I`, and M bases that fall
  in an intron are soft-clipped (so such reads are no longer lost). Works on
  both strands; every output is validated (query length == M+I+S, in-bounds).
  Reads containing a reference deletion (`D`) are still conservatively skipped
  (rare, and not yet remapped). 5 new regression tests.

## [0.0.218] - 2026-08-27

### Fixed (found by a full code re-read)
- **`metagene` options were silently ignored**: `--region`, `--normalize`,
  `--score-transform`, and `--weight-names` were accepted and documented but did
  nothing. They are now implemented: `--region` filters to the chosen gene
  region before profiling; `--normalize` scales each site's weight by its
  transcript length; `--score-transform log2|log10` transforms weight columns;
  `--weight-names` renames the weight columns (so the binned output uses those
  names). Defaults preserve previous behaviour.
- **`mapping` report**: the paired-end "Aligned discordantly 1 time" line was a
  hardcoded `0 (0.02%)` and unpaired/paired line-label semantics were off; the
  report no longer hardcodes a bogus percentage.
- **`mapping` output BAM compression threads** were hardcoded to `4` regardless
  of `--threads`; now uses the requested thread count.
- **`download` `--list` count** of downloaded references counted every `.parquet`
  in the cache dir; now counts only the built-in references present.

## [0.0.217] - 2026-08-27

### Performance
- **`gbam2tbam` is now parallel**: `--threads > 1` uses a `ProcessPoolExecutor`
  like `tbam2gbam` (previously the thread count only affected sorting, so the
  conversion itself ran single-threaded). Also fixed a bug where genome-aligned
  reads were deserialized inside workers against the *transcript* header
  (which broke reference resolution and silently dropped reads); workers now
  deserialize with the correct genome header.
- **`annotate` streams input** instead of `readlines()`-ing the whole file
  (much lower memory on large variant sets), and uses a sorted gene-body
  interval index per chromosome so each site is looked up with `bisect` rather
  than scanning every transcript on the chromosome.

## [0.0.216] - 2026-08-27

### Changed
- **`liftover` is now bidirectional**: `--direction t2g` (transcript -> genome,
  default) or `g2t` (genome -> transcript). `tbam2gbam` and `gbam2tbam` remain
  as clear single-direction convenience commands; the previous duplicate
  `liftover`/`tbam2gbam` help text is now distinct.
- **`annotate` true intronic detection**: a site lying inside a gene's
  transcribed span (between its first and last exon) but not in an exon is now
  classified **Intronic** (assigned to its gene) instead of **Intergenic**.
  Previously intronic detection only fired for positions that were exonic in one
  overlapping transcript yet intronic in another.
- **Comprehensive synthetic-data test suite** (`tests/synthetic_data.py` +
  `tests/test_synthetic_pipeline.py`, 14 integration tests): a small
  hand-computable genome/GTF validates prepare (`--with-txpos`/`--with-codon`),
  annotate (both strands, 5'UTR/CDS/3'UTR/Intronic/Intergenic/Noncoding,
  variant effect with exact codons), motif boundary cases, gbam2tbam, liftover
  direction selection, map_to_local, coordinate and logo - all against exact
  expected values.

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
