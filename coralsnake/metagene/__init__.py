#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Metagene Analysis Package - Using ruranges and polars.
# Migrated from the standalone `metagene` package into the coralsnake CLI
# as the `coralsnake metagene` subcommand. All functions are re-exported here
# so the Python API is identical to the standalone package.

# Core modules with clean API names
from .gtf import prepare_exon_ref, load_gtf
from .io import load_sites, parse_feature_file, load_reference
from .overlap import annotate_with_features, calculate_bin_statistics
from .annotation import map_to_transcripts, normalize_positions, show_summary_stats
from .plotting import plot_profile
from .map_to_local import map_to_local


# Export main functions
__all__ = [
    # Core analysis functions
    "annotate_with_features",
    "calculate_bin_statistics",
    "map_to_transcripts",
    "map_to_local",
    "normalize_positions",
    "show_summary_stats",
    # Data I/O functions
    "load_sites",
    "parse_feature_file",
    "load_reference",
    # GTF parsing functions
    "prepare_exon_ref",
    "load_gtf",
    # Plotting functions
    "plot_profile",
]
