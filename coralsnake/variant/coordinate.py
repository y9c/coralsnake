#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Migrated from the standalone `variant` package (`variant coordinate`).
# urllib3 is replaced with the stdlib urllib; gzip I/O reused from coralsnake.

import os
import sys

from ..utils import get_logger

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


def get_mapper(reference, mapper_type, cache=None):
    """Build a chrom-name mapper from a UCSC chromAlias.txt (cached).

    ``mapper_type`` is ``U2E`` (UCSC→Ensembl) or ``E2U`` (Ensembl→UCSC).
    """
    if reference == "hg38":
        names = ["ucsc", "assembly", "ensembl", "genbank", "refseq"]
    elif reference == "mm39":
        names = ["ucsc", "assembly", "ensembl", "genbank", "refseq"]
    else:
        LOGGER.error(f"Invalid reference: {reference}!")
        sys.exit(1)
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

    def open_file(path):
        return open(path, "r")

    with open_file(reference_path) as mapper_file:
        for line in mapper_file:
            if line.startswith("#"):
                continue
            cols = line.strip("\n").split("\t")
            mapper = dict(zip(names, cols))
            if mapper_type == "U2E":
                chrom_mapper[mapper.get("ucsc", "")] = mapper.get("ensembl", "")
            elif mapper_type == "E2U":
                chrom_mapper[mapper.get("ensembl", "")] = mapper.get("ucsc", "")
    return chrom_mapper


def run_coordinate(
    input_file,
    output_file,
    reference_mapping,
    buildin_mapping,
    columns,
    with_header,
    keep_original,
):
    col_sep = "\t"
    columns_index = [int(x) - 1 for x in str(columns).split(",")]
    if len(columns_index) <= 3:
        columns_index_mapper = dict(zip(["chrom", "pos", "strand"], columns_index))
    else:
        LOGGER.error("Invalid number of columns!")
        sys.exit(1)

    chrom_col = columns_index_mapper.get("chrom")

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
        if buildin_mapping == "U2E":
            chrom_mapper = dict(
                [("chr" + str(i), str(i)) for i in range(1, 100)]
                + [("chrX", "X"), ("chrY", "Y"), ("chrM", "MT")]
            )
        elif buildin_mapping == "E2U":
            chrom_mapper = dict(
                [(str(i), "chr" + str(i)) for i in range(1, 100)]
                + [("X", "chrX"), ("Y", "chrY"), ("MT", "chrM")]
            )
        elif buildin_mapping in ["U2E-hg38", "E2U-hg38", "U2E-mm39", "E2U-mm39"]:
            chrom_mapper = get_mapper(
                buildin_mapping.split("-")[1], buildin_mapping.split("-")[0]
            )
        else:
            LOGGER.error("Invalid buildin_mapping!")
            sys.exit(1)
    else:
        LOGGER.warning("No mapping provided!")
        chrom_mapper = {}

    import gzip

    def open_in(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

    def open_out(path):
        return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

    with open_in(input_file) as input_handle, open_out(output_file) as output_handle:

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

        # read first line and check column number
        input_cols = input_handle.readline().strip("\n").split(col_sep)
        if max(columns_index_mapper.values()) > len(input_cols) - 1:
            LOGGER.error(f"Input file only have {len(input_cols)} columns!")
            sys.exit(1)
        input_handle.seek(0)

        if with_header:
            header_line = input_handle.readline()
            if keep_original:
                header_line = header_line.strip("\n") + col_sep + "RenamedChrom" + "\n"
            output_handle.write(header_line)
        for line in input_handle:
            parse_line(line.strip("\n").split(col_sep))
