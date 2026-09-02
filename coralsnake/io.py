#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Data I/O Module - Handles input/output operations

import polars as pl
import os
from pathlib import Path
from .utils import get_cache_dir
from .config import BUILTIN_REFERENCES


def load_sites(
    input_file_name: str,
    with_header: bool = False,
    meta_col_index: list[int] | None = None,
    separator: str = "\t",
) -> pl.DataFrame:
    """
    Load genomic sites from a file using Polars only.
    Returns:
        Polars DataFrame with processed site information (all input columns
        plus the aliased Chromosome/Start/End/Strand columns)
    """
    df = pl.scan_csv(input_file_name, separator=separator, has_header=with_header)
    colnames = list(df.collect_schema())

    if meta_col_index is None:
        raise ValueError("meta_col_index must be provided")
    if any(i < 0 or i >= len(colnames) for i in meta_col_index):
        raise ValueError(
            f"meta columns reference column {max(meta_col_index) + 1} but the input "
            f"has only {len(colnames)} column(s); pass -m with valid 1-based indices "
            f"(e.g. '1,2,3' for Chrom,Site,Strand or '1,2,3,6' for Chrom,Start,End,Strand)"
        )

    meta_col_names = [colnames[i] for i in meta_col_index]
    # Rename any header named exactly Chromosome/Start/End/Strand so the alias
    # assignment below cannot collide with it.
    newnames = [
        "_original_" + col if col in ["Chromosome", "Start", "End", "Strand"] else col
        for col in colnames
    ]
    # The meta columns' names AFTER the rename (schema_overrides are keyed by
    # the pre-rename names, the with_columns aliases by the post-rename names).
    meta_col_names_renamed = [newnames[i] for i in meta_col_index]

    df = pl.scan_csv(
        input_file_name,
        separator=separator,
        has_header=with_header,
        new_columns=newnames,
        schema_overrides={meta_col_names[0]: pl.Utf8, meta_col_names[-1]: pl.Utf8},
    )

    if len(meta_col_names_renamed) == 4:
        df = df.with_columns(
            pl.col(meta_col_names_renamed[0]).alias("Chromosome"),
            pl.col(meta_col_names_renamed[1]).alias("Start"),
            pl.col(meta_col_names_renamed[2]).alias("End"),
            pl.col(meta_col_names_renamed[3]).alias("Strand"),
        )
    elif len(meta_col_names_renamed) == 3:
        df = df.with_columns(
            pl.col(meta_col_names_renamed[0]).alias("Chromosome"),
            (pl.col(meta_col_names_renamed[1]) - 1).alias("Start"),
            pl.col(meta_col_names_renamed[1]).alias("End"),
            pl.col(meta_col_names_renamed[2]).alias("Strand"),
        )
    else:
        raise ValueError("meta_col_index must specify either 3 or 4 column indices")

    return df.collect()


def parse_feature_file(feature_file_name: str) -> pl.DataFrame:
    """
    Parse a reference file (Parquet) using Polars.

    The reference schema is validated so a corrupt or incompatible cached
    file fails fast with an actionable message instead of a cryptic
    KeyError deep in the analysis.

    Returns:
        Polars DataFrame with processed feature information
    """
    df = pl.read_parquet(feature_file_name)
    required = (
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "gene_id",
        "transcript_id",
        "transcript_length",
        "Start_exon",
        "End_exon",
    )
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Reference file {feature_file_name} is missing required columns "
            f"{missing}; the file may be corrupt or from an incompatible "
            "version. Delete it and re-download it "
            "(coralsnake reference download <ref>)."
        )
    return df


def load_reference(species: str | None = None) -> pl.DataFrame | dict:
    """
    Load built-in reference annotations for common species using Polars only.

    Args:
        species: Species name to load, or None to get available species

    Returns:
        Polars DataFrame with feature annotations, or dict of available species if species=None
    """
    if species is None:
        available = {}
        cache_dir = get_cache_dir()
        for species_name, info in BUILTIN_REFERENCES.items():
            cache_path = cache_dir / Path(info["parquet_file"]).name
            if cache_path.exists():
                file_size_mb = os.path.getsize(cache_path) / (1024 * 1024)
                available[species_name] = {
                    "file": info["parquet_file"],
                    "source": info["source_file"],
                    "description": info["description"],
                    "size_mb": round(file_size_mb, 2),
                    "location": "cache",
                }
        return available

    if species not in BUILTIN_REFERENCES:
        available_species = list(BUILTIN_REFERENCES.keys())
        raise ValueError(
            f"Species '{species}' not available. "
            f"Available species: {available_species}\n"
            f"Use load_reference() without arguments to see available species."
        )

    info = BUILTIN_REFERENCES[species]
    cache_dir = get_cache_dir()
    cache_path = cache_dir / Path(info["parquet_file"]).name

    if cache_path.exists():
        return parse_feature_file(str(cache_path))

    # File doesn't exist - try non-interactive auto-download in CI/pytest/non-tty
    try:
        import click
        from .download import download_references
        import sys

        click.echo(f"Reference '{species}' not found locally.")
        click.echo(f"Description: {info['description']}")
        # If no tty (e.g., under pytest) or METAGENE_AUTO_DOWNLOAD=1, auto-download silently
        auto_env = os.environ.get("METAGENE_AUTO_DOWNLOAD", "").lower() in {
            "1",
            "true",
            "yes",
        }
        if (sys.stdin is None or not sys.stdin.isatty()) or auto_env:
            try:
                download_references(species, silent=True)
                return parse_feature_file(str(cache_path))
            except Exception as e:
                raise RuntimeError(f"Failed to download {species}: {e}")

        if click.confirm("Would you like to download it now?", default=True):
            try:
                download_references(species, silent=True)
                return parse_feature_file(str(cache_path))
            except Exception as e:
                raise RuntimeError(f"Failed to download {species}: {e}")
        else:
            raise ValueError(
                f"Reference '{species}' is required but not available locally."
            )
    except ImportError:
        # click not available, just raise an error
        raise ValueError(
            f"Reference '{species}' not found locally and cannot prompt for download."
        )

    raise ValueError(f"Reference '{species}' not found locally.")
