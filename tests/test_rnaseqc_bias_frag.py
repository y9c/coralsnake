#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the GTF-derived 3' bias and fragment-size metrics."""

import pysam
import pytest

from coralsnake.rnaseqc import run_rnaseqc

# a single long exon (>=600 bp) so 3' bias is computable
GTF = """chr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tgene_id \"G3\"; gene_name \"GENE3\";
chr1\tsrc\texon\t1000\t2000\t.\t+\t.\tgene_id \"G3\"; transcript_id \"T3\"; gene_name \"GENE3\";
"""


@pytest.fixture
def paired_inputs(tmp_path):
    gtf = tmp_path / "anno.gtf"
    gtf.write_text(GTF)

    bam = tmp_path / "pe.bam"
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    )

    def _mk(qname, start, cigar, flag, tlen, mstart):
        a = pysam.AlignedSegment(header)
        a.query_name = qname
        a.flag = flag
        a.reference_name = "chr1"
        a.reference_start = start
        a.cigarstring = cigar
        ql = sum(n for op, n in a.cigartuples if op in (0, 1, 4, 7, 8))
        a.query_sequence = "A" * ql
        a.query_qualities = pysam.qualitystring_to_array("I" * ql)
        a.mapping_quality = 60
        a.set_tag("NM", 0)
        a.next_reference_start = mstart
        a.template_length = tlen
        return a

    reads = []
    # 3 pairs, both mates fully inside exon1 [999,2000) 0-based.
    # mate1 forward covers the 5' window, mate2 reverse covers the 3' window.
    for i in range(3):
        q = f"p{i}"
        m1 = 1149 + i
        m2 = 1749 + i
        # mate1: first-of-pair, forward, proper, mate-reverse (flag 99)
        reads.append(_mk(q, m1, "100M", 99, tlen=700, mstart=m2))
        # mate2: second-of-pair, reverse, proper (flag 147)
        reads.append(_mk(q, m2, "100M", 147, tlen=-700, mstart=m1))
    reads.sort(key=lambda r: (r.reference_start, r.query_name))
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)

    return str(gtf), str(bam)


def test_three_prime_bias_and_fragment_size(paired_inputs, tmp_path):
    gtf, bam = paired_inputs
    outdir = tmp_path / "out"
    m = run_rnaseqc(
        bam, gtf, str(outdir), sample="pe", mapping_quality=0,
        detection_threshold=5, bias_offset=150, bias_window=100,
        bias_gene_length=600,
    )

    # 3 fragment sizes of 700 (both mates of each pair land on the same exon)
    assert m["Average Fragment Length"] == pytest.approx(700)
    assert m["Fragment Length Median"] == pytest.approx(700)

    # 5' and 3' windows are both fully covered -> 3' bias = 0.5
    assert m["Genes used in 3' bias"] == 1
    assert m["Mean 3' bias"] == pytest.approx(0.5)
    assert m["Median 3' bias"] == pytest.approx(0.5)

    # fragment-size histogram is written
    with open(outdir / "pe.fragmentSizes.txt") as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "Fragment Size\tCount"
    assert "700\t3" in lines

    # gene counts: 6 reads (unique) -> G3 detected at threshold 5
    assert m["Genes Detected"] == 1


def test_no_fragment_hist_when_no_pairs(tmp_path):
    # scaffold: single-end reads, but gene too short for 3' bias -> no bias gen
    gtf = tmp_path / "a.gtf"
    gtf.write_text(GTF)
    bam = tmp_path / "a.bam"
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    )
    a = pysam.AlignedSegment(header)
    a.query_name = "r"
    a.flag = 0
    a.reference_name = "chr1"
    a.reference_start = 1350
    a.cigarstring = "50M"
    a.query_sequence = "A" * 50
    a.query_qualities = pysam.qualitystring_to_array("I" * 50)
    a.mapping_quality = 60
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        out.write(a)

    outdir = tmp_path / "o"
    m = run_rnaseqc(str(bam), str(gtf), str(outdir), unpaired=True)
    # short gene (< 600) -> no bias; no pairs -> no fragment metrics
    assert "Genes used in 3' bias" not in m
    assert "Average Fragment Length" not in m
