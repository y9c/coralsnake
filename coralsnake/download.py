#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Reference file download functionality

import os
import shutil
import urllib.error
import urllib.request
from pathlib import Path

from rich.console import Console

from .utils import setup_logger, get_cache_dir, get_file_size, format_file_size
from .config import (
    BUILTIN_REFERENCES,
    GENOME_SIZES_MB,
    GENOME_URLS,
    GITHUB_DOWNLOAD_BASE,
    REFERENCE_GROUPS,
    REFERENCE_SIZES_KB,
)

# Set up logger
logger = setup_logger(__name__)

# Rich console for user-facing output (stderr so piped stdout stays clean).
console = Console(stderr=True)


# Emojis for different statuses
class Emojis:
    CHECK = "✅"
    CROSS = "❌"
    INFO = "ℹ️"
    DOWNLOAD = "⬇️"
    FOLDER = "📁"
    DNA = "🧬"
    MOUSE = "🐭"
    BUG = "🦠"


def list_references(console) -> None:
    """List available built-in references with Rich formatting and emojis."""
    # Get cache directory to check downloaded files
    cache_dir = get_cache_dir()

    # Get list of actually cached files
    cached_files = (
        {f.name for f in cache_dir.glob("*.parquet")} if cache_dir.exists() else set()
    )

    console.print(
        f"\n[bold cyan]{Emojis.DNA} Available Built-in References:[/bold cyan]"
    )
    console.print(f"[dim]{Emojis.FOLDER} Cache directory: {cache_dir}[/dim]\n")

    # Group by species with emojis
    species_groups = {}
    species_emojis = {
        "Human": "👤",
        "Mouse": "🐭",
        "Arabidopsis": "🌱",
        "Rice": "🌾",
        "Zebrafish": "🐟",
        "Fruit Fly": "🪰",
        "Worm": "🪱",
        "Yeast": "🍄",
        "Other": "🧬",
    }

    for ref_name, ref_info in BUILTIN_REFERENCES.items():
        # Check if reference is downloaded by looking for the actual file in cache
        expected_file = f"{ref_name}.parquet"
        status_emoji = Emojis.CHECK if expected_file in cached_files else Emojis.CROSS

        # Classify by species with better logic
        if "Human" in ref_info["description"] or any(
            x in ref_name for x in ["GRCh", "hg"]
        ):
            species = "Human"
        elif "Mouse" in ref_info["description"] or any(
            x in ref_name for x in ["GRCm", "mm", "NCBIM"]
        ):
            species = "Mouse"
        elif "Arabidopsis" in ref_info["description"] or "TAIR" in ref_name:
            species = "Arabidopsis"
        elif "Rice" in ref_info["description"] or "IRGSP" in ref_name:
            species = "Rice"
        elif "Zebrafish" in ref_info["description"] or "GRCz" in ref_name:
            species = "Zebrafish"
        elif (
            any(x in ref_info["description"] for x in ["melanogaster", "Drosophila"])
            or "dm" in ref_name
            or "BDGP" in ref_name
        ):
            species = "Fruit Fly"
        elif (
            any(x in ref_info["description"] for x in ["elegans", "C. elegans"])
            or "ce" in ref_name
            or "WBcel" in ref_name
        ):
            species = "Worm"
        elif any(
            x in ref_info["description"]
            for x in ["cerevisiae", "S. cerevisiae", "pombe", "S. pombe"]
        ) or any(x in ref_name for x in ["sacCer", "R64", "ASM294"]):
            species = "Yeast"
        else:
            species = "Other"

        if species not in species_groups:
            species_groups[species] = []
        species_groups[species].append(
            (
                ref_name,
                ref_info["description"] + _size_text(ref_name, cache_dir, cached_files),
                status_emoji,
            )
        )

    # Calculate the padding for alignment
    max_name_length = max(len(name) for name in BUILTIN_REFERENCES.keys())

    # Print organized by species with emojis
    for species in [
        "Human",
        "Mouse",
        "Zebrafish",
        "Fruit Fly",
        "Worm",
        "Yeast",
        "Arabidopsis",
        "Rice",
        "Other",
    ]:
        if species in species_groups:
            emoji = species_emojis.get(species, "🧬")
            console.print(f"[bold yellow]{emoji} {species}:[/bold yellow]")
            for ref_name, description, status in sorted(species_groups[species]):
                status_color = "green" if status == Emojis.CHECK else "red"
                console.print(
                    f"  [{status_color}]{status}[/{status_color}] [green]{ref_name.ljust(max_name_length + 1)}[/green] - {description}"
                )
            console.print()

    # Count downloaded and total references
    total_refs = len(BUILTIN_REFERENCES)
    downloaded_refs = sum(
        f"{name}.parquet" in cached_files for name in BUILTIN_REFERENCES
    )

    console.print(
        f"[dim]{Emojis.INFO} Total: {total_refs} references available ({downloaded_refs} downloaded)[/dim]"
    )
    console.print(
        f"[dim]{Emojis.DOWNLOAD} Use --download <reference> to download a specific reference[/dim]"
    )


def _resolve_references(request: str) -> list[str]:
    """Resolve a download request: a reference name, a group, or 'all'."""
    key = request.lower()
    if key == "all":
        return list(BUILTIN_REFERENCES.keys())
    if key in REFERENCE_GROUPS:
        return list(REFERENCE_GROUPS[key])
    if request in BUILTIN_REFERENCES:
        return [request]
    raise ValueError(
        f"Reference '{request}' not available. "
        f"Available references: {', '.join(BUILTIN_REFERENCES)}; "
        f"groups: {', '.join(REFERENCE_GROUPS)} or 'all'.\n"
        "Use `coralsnake reference list` to see all options."
    )


def _download_with_progress(url: str, tmp_path: Path, silent: bool = False) -> None:
    """Stream ``url`` to ``tmp_path`` with a progress display (unless silent)."""
    with urllib.request.urlopen(url) as response:
        if silent:
            with open(tmp_path, "wb") as f:
                shutil.copyfileobj(response, f)
            return
        from rich.progress import Progress

        total = int(response.headers.get("Content-Length") or 0)
        with Progress(transient=True) as progress:
            task = progress.add_task("download", total=total or None)
            with open(tmp_path, "wb") as f:
                while True:
                    chunk = response.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    progress.update(task, advance=len(chunk))


def _size_text(ref_name: str, cache_dir: Path, cached_files: set) -> str:
    """Size hint for `--list`: the cached size, or the known download size."""
    fname = f"{ref_name}.parquet"
    if fname in cached_files:
        return f" ({format_file_size(os.path.getsize(cache_dir / fname))})"
    kb = REFERENCE_SIZES_KB.get(ref_name)
    if kb:
        return f" (~{format_file_size(kb * 1024)})"
    return ""


def download_references(
    species: str, silent: bool = False, with_genome: bool = False
) -> None:
    """
    Download reference file(s) for the specified species.

    Args:
        species: Species name or 'all' to download all references
        with_genome: Also fetch the linked genome FASTA for each reference
            (large - see config.GENOME_URLS for sources and sizes).
    """
    # Get cache directory
    cache_dir = get_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Determine which references to download: a reference name, a named
    # group (human/mouse), or "all".
    species_to_download = _resolve_references(species)

    if not silent:
        logger.info(f"Downloading {len(species_to_download)} reference file(s)...")
    try:
        for ref in species_to_download:
            info = BUILTIN_REFERENCES[ref]
            target_path = cache_dir / Path(info["parquet_file"]).name

            if target_path.exists():
                if not silent:
                    logger.info(f"Reference '{ref}' already exists in cache, skipping.")
                continue

            # Construct the download URL using the base URL from config
            download_url = f"{GITHUB_DOWNLOAD_BASE}/{Path(info['parquet_file']).name}"

            if not silent:
                console.print(f"\n{Emojis.DNA} [bold cyan]Reference: {ref}[/bold cyan]")
                console.print(
                    f"{Emojis.INFO} [yellow]Description:[/yellow] {info['description']}"
                )

            # Download to a temp path with progress, then atomically rename,
            # so an interrupted transfer can never leave a corrupt file in
            # the cache.
            logger.info(f"Downloading {Path(info['parquet_file']).name}...")
            tmp_path = target_path.with_name(target_path.name + ".tmp")
            try:
                _download_with_progress(download_url, tmp_path, silent=silent)
            except Exception:
                if tmp_path.exists():
                    tmp_path.unlink()
                raise
            os.replace(tmp_path, target_path)
            if not silent:
                size = format_file_size(get_file_size(target_path))
                logger.info(f"{Emojis.CHECK} Successfully downloaded {ref} ({size})")

    except urllib.error.URLError as e:
        raise RuntimeError(
            f"{Emojis.CROSS} Error downloading reference files: {str(e)}"
        )
    except Exception as e:
        raise RuntimeError(f"{Emojis.CROSS} Error processing reference files: {str(e)}")

    if with_genome:
        if not silent:
            logger.info(
                "Also downloading the linked genome FASTA(s) (large files; "
                "see the printed sizes)..."
            )
        for ref in species_to_download:
            download_genome(ref, silent=silent)

    logger.info(f"{Emojis.CHECK} Reference files download completed!")


# ---------------------------------------------------------------------------
# Genome FASTAs (linked, fetched on demand)
#
# Genome sequences are far too large to ship in the fixed `data` release
# (human genomes are ~3 GB uncompressed). Each reference records a stable
# upstream URL (config.GENOME_URLS). `download_genome` streams the .fa.gz
# into the cache, decompresses it, builds a .fa.fai index, and cross-checks
# the headers against the reference annotation's contig names.
# ---------------------------------------------------------------------------


def genome_cache_dir() -> Path:
    """Directory for cached genome FASTAs (kept separate from the parquets)."""
    d = get_cache_dir() / "genomes"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _genome_paths(ref_name: str) -> tuple[Path, Path, Path]:
    """(compressed, decompressed, index) cache paths for a reference genome."""
    gz = genome_cache_dir() / f"{ref_name}.fa.gz"
    fa = genome_cache_dir() / f"{ref_name}.fa"
    return gz, fa, fa.with_name(fa.name + ".fai")


def _gunzip_to(src: Path, dst: Path) -> None:
    """Decompress ``src`` (.gz) to ``dst`` via a temp file (atomic rename)."""
    import gzip

    tmp = dst.with_name(dst.name + ".tmp")
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    os.replace(tmp, dst)


def _check_genome_headers(ref_name: str, fa_path: Path) -> set[str] | None:
    """Cross-check genome FASTA headers against the reference's contigs.

    Raises RuntimeError if there is no overlap at all (wrong genome / naming
    scheme); returns the set of contigs without sequence (None if the
    reference parquet is not cached and the check was skipped).
    """
    import polars as pl

    cache_dir = get_cache_dir()
    info = BUILTIN_REFERENCES[ref_name]
    pq = cache_dir / Path(info["parquet_file"]).name
    if not pq.exists():
        logger.info(f"{ref_name}: reference parquet not cached; skipping header check.")
        return None

    chroms = set(
        pl.scan_parquet(pq)
        .select("Chromosome")
        .unique()
        .collect()["Chromosome"]
        .to_list()
    )
    headers = set()
    with open(fa_path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                headers.add(line[1:].split(None, 1)[0])

    if not chroms & headers:
        raise RuntimeError(
            f"Genome FASTA for '{ref_name}' shares no contig names with the "
            f"reference annotation ({len(chroms)} contigs, e.g. "
            f"{', '.join(sorted(chroms)[:3])}); the linked FASTA "
            f"({GENOME_URLS[ref_name]}) does not match this annotation."
        )
    missing = chroms - headers
    if missing:
        sample = ", ".join(sorted(missing)[:5])
        more = f" (+{len(missing) - 5} more)" if len(missing) > 5 else ""
        logger.warning(
            f"{ref_name}: {len(missing)} annotated contig(s) have no sequence "
            f"in the linked genome FASTA (e.g. {sample}{more})"
        )
    return missing


def download_genome(ref_name: str, silent: bool = False) -> Path:
    """Fetch (and cache) the linked genome FASTA for a built-in reference.

    Returns the path to the decompressed FASTA (the ``.fa.fai`` index is
    built alongside). Each step is idempotent: download, decompress and
    indexing only run when their artifact is missing.
    """
    info = BUILTIN_REFERENCES.get(ref_name)
    if info is None:
        raise ValueError(
            f"Unknown reference '{ref_name}'. "
            f"Available: {', '.join(BUILTIN_REFERENCES)}."
        )
    if ref_name not in GENOME_URLS:
        raise ValueError(f"No genome FASTA link recorded for '{ref_name}'.")

    gz, fa, fai = _genome_paths(ref_name)

    if not silent:
        console.print(f"\n{Emojis.DNA} [bold cyan]Genome: {ref_name}[/bold cyan]")
        console.print(
            f"{Emojis.INFO} [yellow]Description:[/yellow] {info['description']}"
        )
        mb = GENOME_SIZES_MB.get(ref_name)
        if mb:
            console.print(f"{Emojis.DOWNLOAD} [yellow]~{mb} MB compressed[/yellow]")
        console.print(f"{Emojis.INFO} [dim]{GENOME_URLS[ref_name]}[/dim]")

    if not gz.exists():
        logger.info(f"Downloading genome FASTA for '{ref_name}'...")
        tmp = gz.with_name(gz.name + ".tmp")
        try:
            _download_with_progress(GENOME_URLS[ref_name], tmp, silent=silent)
        except Exception:
            if tmp.exists():
                tmp.unlink()
            raise
        os.replace(tmp, gz)

    if not fa.exists():
        logger.info(f"Decompressing {gz.name} (this can take a minute)...")
        _gunzip_to(gz, fa)

    if not fai.exists():
        import pysam

        logger.info("Building FASTA index (.fai)...")
        pysam.faidx(str(fa))

    _check_genome_headers(ref_name, fa)

    if not silent:
        logger.info(f"{Emojis.CHECK} Genome '{ref_name}' ready: {fa}")
    return fa


# ---------------------------------------------------------------------------
# Reference name -> cached artifact (parquet / derived table / derived GTF)
#
# Lets the analysis tools take a reference NAME instead of a file path: the
# name resolves to the shared cached object, downloading the parquet (or
# deriving the table/GTF) on first use. Derived files sit next to the
# parquet in the cache and are re-exported when the parquet is newer.
# ---------------------------------------------------------------------------


def ensure_reference(ref_name: str, silent: bool = False) -> Path:
    """Return the cached parquet path for a built-in reference, downloading
    it first if it is missing."""
    info = BUILTIN_REFERENCES.get(ref_name)
    if info is None:
        raise ValueError(
            f"Unknown reference '{ref_name}'. "
            f"Available: {', '.join(BUILTIN_REFERENCES)}."
        )
    cache_dir = get_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / Path(info["parquet_file"]).name
    if not target.exists():
        download_references(ref_name, silent=silent)
    if not target.exists():
        raise RuntimeError(
            f"Reference parquet {target.name} is still missing after the "
            "download attempt."
        )
    return target


def _derived_artifact(ref_name: str, suffix: str) -> Path:
    return get_cache_dir() / f"{ref_name}.{suffix}"


def _ensure_derived(ref_name: str, suffix: str, export, silent: bool) -> Path:
    from .io import parse_feature_file

    pq = ensure_reference(ref_name, silent=silent)
    out = _derived_artifact(ref_name, suffix)
    if out.exists() and out.stat().st_mtime >= pq.stat().st_mtime:
        return out
    df = parse_feature_file(str(pq))
    tmp = out.with_name(out.name + ".tmp")
    export(df, str(tmp))
    os.replace(tmp, out)
    return out


def ensure_reference_table(ref_name: str, silent: bool = False) -> Path:
    """Cached ``prepare``-style annotation table (TSV) derived from the
    reference parquet; feeds ``liftover -a`` and ``annotate --annotation``."""
    from .ref_export import export_table

    return _ensure_derived(ref_name, "table.tsv", export_table, silent)


def ensure_reference_gtf(ref_name: str, silent: bool = False) -> Path:
    """Cached GTF derived from the reference parquet; feeds
    ``annotate --reference-gtf``."""
    from .ref_export import export_gtf

    return _ensure_derived(ref_name, "gtf", export_gtf, silent)
