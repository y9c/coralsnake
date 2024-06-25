#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-05-25 16:14


import os
import subprocess

from .conversion import convert_file


def run_shell_command(command):
    """Utility function to run a shell command."""
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error executing command: {command}\nError: {stderr.decode()}")
    else:
        print(f"Command executed successfully: {command}\nOutput: {stdout.decode()}")
    return process.returncode


def combine_genes_fa(config, INTERNALDIR):
    input_files = [
        os.path.expanduser(f)
        for f in config.get("spike", [])
        + [INTERNALDIR / "reference_file/transcript.fa"]
    ]
    output_fa = INTERNALDIR / "reference_file/genes.fa"
    output_fai = INTERNALDIR / "reference_file/genes.fa.fai"

    cat_command = f"cat {' '.join(input_files)} > {output_fa}"
    run_shell_command(cat_command)

    faidx_command = f"samtools faidx {output_fa} --fai-idx {output_fai}"
    run_shell_command(faidx_command)


def make_ref_conversion(INTERNALDIR):
    input_fa = INTERNALDIR / "reference_file/genes.fa"
    output_mk = INTERNALDIR / "reference_file/ref.MK.fa"
    output_km = INTERNALDIR / "reference_file/ref.KM.fa"
    convert_file(input_fa, output_mk, output_km)


def build_bwa_index(INTERNALDIR, type_):
    input_fa = INTERNALDIR / f"reference_file/ref.{type_}.fa"
    params = str(INTERNALDIR / f"reference_file/ref.{type_}")

    command = f"~/tools/bwa-mem2/bwa-mem2 index -p {params} {input_fa}"
    run_shell_command(command)


# The following functions (map_to_fwd_PE, map_to_rev_PE, map_to_fwd_SE, map_to_rev_SE, combine_map_join, replace_reads, combine_map_sortname, replace_reads_calmd, add_tag, filter_and_sort_bam, select_unmapped_reads_PE, select_unmapped_reads_SE) need to be defined similarly as the above functions.

# Example usage:
# config = {...}
# INTERNALDIR = "path_to_internaldir"
# SAMPLE2DATA = {...}
# SAMPLE2LIB = {...}
# TEMPDIR = "path_to_tempdir"
# threads = 16
# combine_genes_fa(config, INTERNALDIR)
# make_ref_conversion(INTERNALDIR)
# build_bwa_index(INTERNALDIR, "MK")
# cutadapt_SE(SAMPLE2DATA, TEMPDIR, INTERNALDIR, "sample1", "rn1", threads)
# cutadapt_PE(SAMPLE2DATA, TEMPDIR, INTERNALDIR, "sample1", "rn1", threads)
# report_cutadapt(SAMPLE2DATA)
# qc_trimmed(SAMPLE2DATA, TEMPDIR, "sample1", "rn1", "R1")
# report_qc_trimmed(SAMPLE2DATA)
# convert_reads(SAMPLE2DATA, TEMPDIR, "sample1", "rn1", "R1")
