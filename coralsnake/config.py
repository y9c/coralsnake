#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Configuration data for Metagene

from typing import TypedDict


class ReferenceInfo(TypedDict):
    parquet_file: str
    source_file: str
    source_url: str
    description: str


# GitHub repository information.
#
# The reference parquets are served from this repository's `data` release
# (prerelease). That release is FIXED: its assets are immutable, so a data
# update is published under a new tag (e.g. `data-v2`) together with a bump
# of GITHUB_DOWNLOAD_BASE here. Files are (re)built by scripts/build_references.py.
GITHUB_REPO = "y9c/coralsnake"
GITHUB_API_URL = f"https://api.github.com/repos/{GITHUB_REPO}/releases/latest"
GITHUB_DOWNLOAD_BASE = f"https://github.com/{GITHUB_REPO}/releases/download/data"

# Download groups for `coralsnake reference download <group>`
# ("all" is handled implicitly: every reference in BUILTIN_REFERENCES.)
REFERENCE_GROUPS: dict[str, list[str]] = {
    "human": ["GRCh38", "hg38", "GRCh37", "hg19"],
    "mouse": ["GRCm39", "mm39", "GRCm38", "mm10", "mm9", "NCBIM37"],
}

# Approximate download size of each hosted parquet (KB). Shown by
# `coralsnake reference list`. Keep in sync when the data release is
# updated (scripts/build_references.py prints the new sizes on every build).
REFERENCE_SIZES_KB: dict[str, int] = {
    "GRCh38": 15411,
    "hg38": 20306,
    "GRCh37": 12966,
    "hg19": 19122,
    "GRCm39": 9188,
    "mm39": 11627,
    "GRCm38": 8545,
    "mm10": 11234,
    "mm9": 6452,
    "NCBIM37": 7105,
    "TAIR10": 3598,
    "IRGSP-1.0": 2454,
    "WBcel235": 3591,
    "ce11": 2807,
    "BDGP6.32": 2222,
    "dm6": 2252,
    "GRCz11": 6300,
    "GRCz10": 6617,
    "bGalGal1": 6039,
    "Glycine_max_v2.1": 5662,
    "R64-1-1": 119,
    "sacCer3": 123,
    "ASM294v2": 232,
}

# Genome FASTA links (NOT shipped in the data release - far too large to
# host; human genomes are ~3 GB uncompressed each). ``reference genome <ref>``
# (or ``reference download <ref> --with-genome``) streams the file into the
# cache on demand, decompresses it, builds a .fa.fai index and cross-checks
# the headers against the reference's contig names. URLs verified
# 2026-09-02 (HTTP 200); headers match each reference parquet's Chromosome
# values. Note the BDGP6.32 DNA is taken from release-109: release-110 moved
# D. melanogaster to assembly BDGP6.46, and the genome sequence of a given
# assembly is identical across releases.
GENOME_URLS: dict[str, str] = {
    "GRCh38": "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    "GRCh37": "https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz",
    "GRCm39": "https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
    "GRCm38": "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz",
    "NCBIM37": "https://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz",
    "TAIR10": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
    "IRGSP-1.0": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz",
    "WBcel235": "https://ftp.ensembl.org/pub/release-110/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz",
    "BDGP6.32": "https://ftp.ensembl.org/pub/release-109/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz",
    "GRCz11": "https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz",
    "GRCz10": "https://ftp.ensembl.org/pub/release-91/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz",
    "bGalGal1": "https://ftp.ensembl.org/pub/release-110/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz",
    "Glycine_max_v2.1": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz",
    "R64-1-1": "https://ftp.ensembl.org/pub/release-57/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.SGD1.01.57.dna.toplevel.fa.gz",
    "ASM294v2": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi/current/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz",
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "hg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    "mm39": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
    "mm10": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
    "mm9": "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz",
    "ce11": "https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz",
    "dm6": "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz",
    "sacCer3": "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
}

# Approximate COMPRESSED download size of each linked genome FASTA (MB),
# shown before a genome download. Measured 2026-09-02 (Content-Length).
GENOME_SIZES_MB: dict[str, int] = {
    "GRCh38": 842,
    "GRCh37": 830,
    "GRCm39": 769,
    "GRCm38": 769,
    "NCBIM37": 729,
    "TAIR10": 35,
    "IRGSP-1.0": 109,
    "WBcel235": 29,
    "BDGP6.32": 41,
    "GRCz11": 391,
    "GRCz10": 391,
    "bGalGal1": 304,
    "Glycine_max_v2.1": 274,
    "R64-1-1": 4,
    "ASM294v2": 4,
    "hg38": 939,
    "hg19": 904,
    "mm39": 831,
    "mm10": 830,
    "mm9": 821,
    "ce11": 30,
    "dm6": 43,
    "sacCer3": 4,
}

# Built-in reference mappings with source information
#
# source_file: local layout (relative to the source dir used by
#   scripts/build_references.py) - "<Species>/raw/<upstream filename>".
# source_url: verified upstream location (checked 2026-09-02). Note the
#   upstream split: animals + S. cerevisiae live on ftp.ensembl.org,
#   plants + S. pombe on ftp.ebi.ac.uk/ensemblgenomes (Ensembl Genomes),
#   and UCSC gene GTFs moved under bigZips/genes/ (the old goldenPath
#   "genes/" dirs and dated combined GTFs are gone).
BUILTIN_REFERENCES: dict[str, ReferenceInfo] = {
    # Human genomes
    "GRCh38": {
        "parquet_file": "GRCh38.parquet",
        "source_file": "Homo_sapiens/raw/Homo_sapiens.GRCh38.110.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz",
        "description": "Human genome GRCh38 (Ensembl release 110)",
    },
    "hg38": {
        "parquet_file": "hg38.parquet",
        "source_file": "Homo_sapiens/raw/hg38.knownGene.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz",
        "description": "Human genome hg38 (UCSC GENCODE/knownGene)",
    },
    "GRCh37": {
        "parquet_file": "GRCh37.parquet",
        "source_file": "Homo_sapiens/raw/Homo_sapiens.GRCh37.75.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz",
        "description": "Human genome GRCh37 (Ensembl release 75)",
    },
    "hg19": {
        "parquet_file": "hg19.parquet",
        "source_file": "Homo_sapiens/raw/hg19.knownGene.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.knownGene.gtf.gz",
        "description": "Human genome hg19 (UCSC GENCODE/knownGene)",
    },
    # Mouse genomes
    "GRCm39": {
        "parquet_file": "GRCm39.parquet",
        "source_file": "Mus_musculus/raw/Mus_musculus.GRCm39.110.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz",
        "description": "Mouse genome GRCm39 (Ensembl release 110)",
    },
    "mm39": {
        "parquet_file": "mm39.parquet",
        "source_file": "Mus_musculus/raw/mm39.knownGene.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.knownGene.gtf.gz",
        "description": "Mouse genome mm39 (UCSC GENCODE/knownGene)",
    },
    "GRCm38": {
        "parquet_file": "GRCm38.parquet",
        "source_file": "Mus_musculus/raw/Mus_musculus.GRCm38.102.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz",
        "description": "Mouse genome GRCm38 (Ensembl release 102)",
    },
    "mm10": {
        "parquet_file": "mm10.parquet",
        "source_file": "Mus_musculus/raw/mm10.knownGene.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.knownGene.gtf.gz",
        "description": "Mouse genome mm10 (UCSC GENCODE/knownGene)",
    },
    "mm9": {
        "parquet_file": "mm9.parquet",
        "source_file": "Mus_musculus/raw/mm9.knownGene.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/genes/mm9.knownGene.gtf.gz",
        "description": "Mouse genome mm9 (UCSC GENCODE/knownGene)",
    },
    "NCBIM37": {
        "parquet_file": "NCBIM37.parquet",
        "source_file": "Mus_musculus/raw/Mus_musculus.NCBIM37.67.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz",
        "description": "Mouse genome NCBIM37 (Ensembl release 67)",
    },
    # Model organisms
    "TAIR10": {
        "parquet_file": "TAIR10.parquet",
        "source_file": "Arabidopsis_thaliana/raw/Arabidopsis_thaliana.TAIR10.63.gtf.gz",
        "source_url": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gtf.gz",
        "description": "Arabidopsis thaliana TAIR10 (EnsemblGenomes r63)",
    },
    "IRGSP-1.0": {
        "parquet_file": "IRGSP-1.0.parquet",
        "source_file": "Oryza_sativa/raw/Oryza_sativa.IRGSP-1.0.63.gtf.gz",
        "source_url": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.63.gtf.gz",
        "description": "Rice IRGSP-1.0 (EnsemblGenomes r63)",
    },
    "WBcel235": {
        "parquet_file": "WBcel235.parquet",
        "source_file": "Caenorhabditis_elegans/raw/Caenorhabditis_elegans.WBcel235.110.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-110/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.110.gtf.gz",
        "description": "C. elegans WBcel235 (Ensembl release 110)",
    },
    "ce11": {
        "parquet_file": "ce11.parquet",
        "source_file": "Caenorhabditis_elegans/raw/ce11.refGene.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.refGene.gtf.gz",
        "description": "C. elegans ce11 (UCSC refGene/WormBase)",
    },
    "BDGP6.32": {
        "parquet_file": "BDGP6.32.parquet",
        "source_file": "Drosophila_melanogaster/raw/Drosophila_melanogaster.BDGP6.32.110.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.110.gtf.gz",
        "description": "D. melanogaster BDGP6.32 (Ensembl release 110)",
    },
    "dm6": {
        "parquet_file": "dm6.parquet",
        "source_file": "Drosophila_melanogaster/raw/dm6.refGene.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.refGene.gtf.gz",
        "description": "D. melanogaster dm6 (UCSC refGene)",
    },
    "GRCz11": {
        "parquet_file": "GRCz11.parquet",
        "source_file": "Danio_rerio/raw/Danio_rerio.GRCz11.110.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-110/gtf/danio_rerio/Danio_rerio.GRCz11.110.gtf.gz",
        "description": "Zebrafish GRCz11 (Ensembl release 110)",
    },
    "GRCz10": {
        "parquet_file": "GRCz10.parquet",
        "source_file": "Danio_rerio/raw/Danio_rerio.GRCz10.91.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-91/gtf/danio_rerio/Danio_rerio.GRCz10.91.gtf.gz",
        "description": "Zebrafish GRCz10 (Ensembl release 91)",
    },
    "bGalGal1": {
        "parquet_file": "bGalGal1.parquet",
        "source_file": "Gallus_gallus/raw/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-110/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf.gz",
        "description": "Chicken bGalGal1.mat.broiler.GRCg7b (Ensembl release 110)",
    },
    "Glycine_max_v2.1": {
        "parquet_file": "Glycine_max_v2.1.parquet",
        "source_file": "Glycine_max/raw/Glycine_max.Glycine_max_v2.1.63.gtf.gz",
        "source_url": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/gtf/glycine_max/Glycine_max.Glycine_max_v2.1.63.gtf.gz",
        "description": "Soybean Glycine max v2.1 (EnsemblGenomes r63)",
    },
    "R64-1-1": {
        "parquet_file": "R64-1-1.parquet",
        "source_file": "Saccharomyces_cerevisiae/raw/Saccharomyces_cerevisiae.SGD1.01.57.gtf.gz",
        "source_url": "https://ftp.ensembl.org/pub/release-57/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.SGD1.01.57.gtf.gz",
        "description": "S. cerevisiae R64-1-1 (Ensembl release 57)",
    },
    "sacCer3": {
        "parquet_file": "sacCer3.parquet",
        "source_file": "Saccharomyces_cerevisiae/raw/sacCer3.ncbiRefSeq.gtf.gz",
        "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ncbiRefSeq.gtf.gz",
        "description": "S. cerevisiae sacCer3 (UCSC ncbiRefSeq/SGD)",
    },
    "ASM294v2": {
        "parquet_file": "ASM294v2.parquet",
        "source_file": "Schizosaccharomyces_pombe/raw/Schizosaccharomyces_pombe.ASM294v2.63.gtf.gz",
        "source_url": "https://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi/current/gtf/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.63.gtf.gz",
        "description": "S. pombe ASM294v2 (EnsemblGenomes r63)",
    },
}
