#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Chromosome-name coordinate mapping, fused into coralsnake from the
# standalone `variant` package (`variant coordinate`). urllib3 is replaced by
# the stdlib urllib; gzip/stdin-stdout I/O uses xopen.

import os

from .utils import get_logger

LOGGER = get_logger(__name__)


def download_file(url, path):
    """Stream ``url`` to ``path`` using stdlib urllib (replaces urllib3)."""
    import urllib.request

    with urllib.request.urlopen(url) as r, open(path, "wb") as f:
        while True:
            chunk = r.read(1 << 16)
            if not chunk:
                break
            f.write(chunk)


# Reference name sets for the built-in chromAlias mappers.
_ALIAS_NAMES = ["ucsc", "assembly", "ensembl", "genbank", "refseq"]


def get_mapper(reference, mapper_type, cache=None):
    """Build a chrom-name mapper from a UCSC chromAlias.txt (cached).

    ``mapper_type`` is ``U2E`` (UCSC→Ensembl) or ``E2U`` (Ensembl→UCSC).
    """
    if reference not in ("hg38", "mm39"):
        raise ValueError(f"Invalid reference: {reference}!")

    chrom_mapper = {}

    if cache is None:
        cache = os.path.expanduser("~/.cache/coralsnake/coordinate")
    os.makedirs(cache, exist_ok=True)

    basename = f"{reference}.chromAlias.txt"
    reference_path = os.path.join(cache, basename)
    if not os.path.exists(reference_path):
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/{}/bigZips/{}".format(
            reference,
            basename if reference != "hg38" else "latest/" + basename,
        )
        download_file(url, reference_path)

    with open(reference_path, "r") as mapper_file:
        for line in mapper_file:
            if line.startswith("#"):
                continue
            cols = line.strip("\n").split("\t")
            mapper = dict(zip(_ALIAS_NAMES, cols))
            if mapper_type == "U2E":
                chrom_mapper[mapper.get("ucsc", "")] = mapper.get("ensembl", "")
            elif mapper_type == "E2U":
                chrom_mapper[mapper.get("ensembl", "")] = mapper.get("ucsc", "")
    return chrom_mapper


def _builtin_mapper(buildin_mapping):
    """Return a built-in chrom-name mapper for a preset string."""
    if buildin_mapping in ("U2E", "E2U"):
        if buildin_mapping == "U2E":
            return dict(
                [("chr" + str(i), str(i)) for i in range(1, 100)]
                + [("chrX", "X"), ("chrY", "Y"), ("chrM", "MT")]
            )
        return dict(
            [(str(i), "chr" + str(i)) for i in range(1, 100)]
            + [("X", "chrX"), ("Y", "chrY"), ("MT", "chrM")]
        )
    if buildin_mapping in ("U2E-hg38", "E2U-hg38", "U2E-mm39", "E2U-mm39"):
        return get_mapper(buildin_mapping.split("-")[1], buildin_mapping.split("-")[0])
    raise ValueError(f"Invalid buildin_mapping: {buildin_mapping}!")


def run_coordinate(
    input_file,
    output_file,
    reference_mapping,
    buildin_mapping,
    columns,
    with_header,
    keep_original,
):
    """Rename the chrom column of every input row."""
    from xopen import xopen

    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    if len(columns_index) > 3:
        raise ValueError("Invalid number of columns!")
    if len(columns_index) < 1:
        raise ValueError("Need at least one column (the chrom column)!")
    chrom_col = columns_index[0]

    if reference_mapping:
        if buildin_mapping:
            LOGGER.warning(
                "Both reference_mapping and buildin_mapping are provided, "
                "reference_mapping will be used!"
            )
        with open(reference_mapping) as mapper_file:
            chrom_mapper = dict(
                (line.strip("\n").split("\t")[:2] for line in mapper_file)
            )
    elif buildin_mapping:
        chrom_mapper = _builtin_mapper(buildin_mapping)
    else:
        LOGGER.warning("No mapping provided!")
        chrom_mapper = {}

    with xopen(input_file, "rt") as input_handle, xopen(output_file, "wt") as output_handle:

        def parse_line(input_cols):
            chrom = input_cols[chrom_col]
            chrom_rename = chrom_mapper.get(chrom, chrom)
            if keep_original:
                output_cols = input_cols + [chrom_rename]
            else:
                output_cols = (
                    input_cols[:chrom_col]
                    + [chrom_rename]
                    + input_cols[chrom_col + 1 :]
                )
            output_handle.write(col_sep.join(output_cols) + "\n")

        # Read the first line to decide header / row layout (no seek, so stdin works).
        first = input_handle.readline()
        if not first:
            return  # empty input
        if max(columns_index) > len(first.rstrip("\n").split(col_sep)) - 1:
            raise ValueError(
                f"Input file only has {len(first.rstrip().split(col_sep))} columns!"
            )

        if with_header:
            header_line = first
            if keep_original:
                header_line = header_line.strip("\n") + col_sep + "RenamedChrom" + "\n"
            output_handle.write(header_line)
        else:
            parse_line(first.rstrip("\n").split(col_sep))
        for line in input_handle:
            parse_line(line.strip("\n").split(col_sep))
