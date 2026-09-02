#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``coralsnake.rnaseqc`` (the ``qc`` RNA-seq metrics).

Builds a tiny synthetic GTF + BAM in ``tmp_path`` and checks the exact QC
metrics (mapping statistics, exon/intronic/intergenic/ambiguous classification,
gene counts / TPM, genes detected).
"""

import pysam
import pytest

from coralsnake.rnaseqc import run_rnaseqc, _estimate_library_complexity

GENE1 = "G1"
GENE2 = "G2"

GTF = """chr1\tsrc\tgene\t101\t400\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";
chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"GENE1\";
chr1\tsrc\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"GENE1\";
chr1\tsrc\tgene\t501\t600\t.\t-\t.\tgene_id \"G2\"; gene_name \"GENE2\";
chr1\tsrc\texon\t501\t600\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\"; gene_name \"GENE2\";
"""


@pytest.fixture
def qc_inputs(tmp_path):
    gtf = tmp_path / "anno.gtf"
    gtf.write_text(GTF)

    bam = tmp_path / "reads.bam"
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    )

    def _mk(qname, start, cigar, flag=0):
        a = pysam.AlignedSegment(header)
        a.query_name = qname
        a.flag = flag
        a.reference_name = "chr1"
        a.reference_start = start
        a.cigarstring = cigar
        # query length = number of bases the CIGAR consumes from the query
        qlen = sum(length for op, length in a.cigartuples if op in (0, 1, 4, 7, 8))
        a.query_sequence = "A" * qlen
        a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
        a.mapping_quality = 60
        a.set_tag("NM", 0)
        return a

    reads = [
        _mk("ex1", 100, "50M"),                 # exonic (G1 exon1)
        _mk("spliced", 100, "25M200N25M"),      # exonic (G1, 2 blocks)
        _mk("ambig", 100, "25M375N25M"),        # ambiguous (G1 exon + G2 exon)
        _mk("ex2", 115, "50M"),                 # exonic (G1 exon1)
        _mk("dup", 130, "30M", flag=0x400),     # exonic + duplicate
        _mk("intr", 220, "20M"),                # intronic (G1 body)
        _mk("inter", 5000, "20M"),              # intergenic
    ]
    reads.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)

    return str(gtf), str(bam)


def _read_metrics(path):
    d = {}
    for line in open(path):
        if "\t" in line:
            k, _, v = line.rstrip("\n").partition("\t")
            d[k] = v
    return d


def test_core_metrics(qc_inputs, tmp_path, capsys):
    gtf, bam = qc_inputs
    outdir = tmp_path / "out"
    m = run_rnaseqc(
        bam, gtf, str(outdir), sample="s", unpaired=True,
        mapping_quality=0, detection_threshold=3,
    )

    assert m["Mapping Rate"] == 1.0
    assert m["Unique Rate of Mapped"] == pytest.approx(6 / 7)
    assert m["Duplicate Rate of Mapped"] == pytest.approx(1 / 7)
    assert m["Exonic Rate"] == pytest.approx(4 / 7)
    assert m["Intronic Rate"] == pytest.approx(1 / 7)
    assert m["Intergenic Rate"] == pytest.approx(1 / 7)
    assert m["Ambiguous Alignment Rate"] == pytest.approx(1 / 7)
    assert m["Intragenic Rate"] == pytest.approx(5 / 7)
    assert m["rRNA Rate"] == 0.0
    assert m["Read Length"] == 50
    assert m["Genes Detected"] == 1
    assert m["Total Reads (Raw)"] == 7

    # output files
    metrics_file = outdir / "s.metrics.tsv"
    assert metrics_file.exists()
    md = _read_metrics(metrics_file)
    assert md["Mapping Rate"] == "1"

    # gene counts (high-quality, unambiguous reads; dup counts but not unique)
    gene_reads = {}
    for line in open(outdir / "s.gene_reads.tsv").read().splitlines()[1:]:
        n, name, c = line.split("\t")
        gene_reads[n] = int(c)
    assert gene_reads["G1"] == 4

    tpm = {}
    for line in open(outdir / "s.gene_tpm.tsv").read().splitlines()[1:]:
        n, name, v = line.split("\t")
        tpm[n] = float(v)
    assert tpm["G1"] == pytest.approx(1e6, rel=1e-3)  # only expressed gene -> 1e6
    assert tpm.get("G2", 0.0) == 0.0  # G2 has no counts -> not emitted


def test_metrics_written_and_rerunnable(qc_inputs, tmp_path):
    gtf, bam = qc_inputs
    outdir = tmp_path / "out2"
    run_rnaseqc(bam, gtf, str(outdir), unpaired=True)
    assert (outdir / "reads.metrics.tsv").exists()
    assert (outdir / "reads.gene_reads.tsv").exists()
    assert (outdir / "reads.gene_tpm.tsv").exists()
    assert (outdir / "reads.exon_reads.tsv").exists()


def test_library_complexity():
    assert _estimate_library_complexity(0, 0) == 0
    assert _estimate_library_complexity(10, 0) == 0
    # more duplication -> LOWER estimated library complexity
    assert _estimate_library_complexity(1000, 500) < _estimate_library_complexity(
        1000, 100
    )
