#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Reference file download functionality

import urllib.request
import urllib.error
from pathlib import Path

from rich.console import Console

from .utils import setup_logger, get_cache_dir, get_file_size, format_file_size
from .config import BUILTIN_REFERENCES, GITHUB_DOWNLOAD_BASE

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
            (ref_name, ref_info["description"], status_emoji)
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
    downloaded_refs = len(cached_files)

    console.print(
        f"[dim]{Emojis.INFO} Total: {total_refs} references available ({downloaded_refs} downloaded)[/dim]"
    )
    console.print(
        f"[dim]{Emojis.DOWNLOAD} Use --download <reference> to download a specific reference[/dim]"
    )


def download_references(species: str, silent: bool = False) -> None:
    """
    Download reference file(s) for the specified species.

    Args:
        species: Species name or 'all' to download all references
    """
    if species not in BUILTIN_REFERENCES and species.lower() != "all":
        raise ValueError(f"Unknown species: {species}")

    # Get cache directory
    cache_dir = get_cache_dir()

    # Determine which species to download
    species_to_download = (
        BUILTIN_REFERENCES.keys() if species.lower() == "all" else [species]
    )

    # Download and process the release
    logger.info("Downloading reference files...")
    try:
        for species in species_to_download:
            info = BUILTIN_REFERENCES[species]
            target_path = cache_dir / Path(info["parquet_file"]).name

            if target_path.exists():
                logger.info(f"Reference for {species} already exists in cache.")
                continue

            # Construct the download URL using the base URL from config
            download_url = f"{GITHUB_DOWNLOAD_BASE}/{Path(info['parquet_file']).name}"

            # Prompt user to download if file doesn't exist
            if not target_path.exists():
                if silent:
                    # Silent mode - just download without prompting
                    pass
                else:
                    # Interactive mode - prompt user via rich console
                    console.print(
                        f"\n{Emojis.DNA} [bold cyan]Reference: {species}[/bold cyan]"
                    )
                    console.print(
                        f"{Emojis.INFO} [yellow]Description:[/yellow] {info['description']}"
                    )
                    console.print(
                        f"{Emojis.DOWNLOAD} [yellow]Source:[/yellow] {info['source_file']}"
                    )
                    response = console.input(
                        f"{Emojis.DOWNLOAD} Reference file not found in cache. "
                        "Download it? (Y/n): "
                    )
                    if response.strip().lower() in ["n", "no"]:
                        logger.info(f"Skipping download for {species}.")
                        continue

            # Download the file
            logger.info(f"Downloading {Path(info['parquet_file']).name}...")
            urllib.request.urlretrieve(download_url, target_path)
            size = format_file_size(get_file_size(target_path))
            logger.info(
                f"{Emojis.CHECK} Successfully downloaded {Path(info['parquet_file']).name} ({size})"
            )

    except urllib.error.URLError as e:
        raise RuntimeError(
            f"{Emojis.CROSS} Error downloading reference files: {str(e)}"
        )
    except Exception as e:
        raise RuntimeError(f"{Emojis.CROSS} Error processing reference files: {str(e)}")

    logger.info(f"{Emojis.CHECK} Reference files download completed!")
