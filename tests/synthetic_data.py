#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Synthetic test data: a small, hand-computable genome + GTF.

Chromosome ``chr1`` is a short, deterministic sequence so every coordinate,
codon and motif can be predicted exactly.

Gene layout (0-based half-open [start, end)):
  g1   '+' protein_coding   exon1 [10,30)  exon2 [50,70)  (t1)
       start_codon at [20,23), stop_codon at [62,65)
          -> CDS = t[10, 52+2) in transcript coords (see fixture)
  g2   '-' protein_coding   exon   [120,150)  (t2)
       start_codon at [137,140), stop_codon at [128,131)  (5' = right end)
  g3   '+' lincRNA          exon   [200,215)  (t3, non-coding -> no codons)
"""

CHR1 = ("ATGCTAGCTAG" * 30)  # 300 bp, deterministic

def _spans_row(gene, tx, chrom, strand, spans, **extra):
    cols = {"gene_id": gene, "transcript_id": tx, "chrom": chrom,
            "strand": strand, "spans": spans}
    cols.update(extra)
    return "\t".join(str(cols[k]) for k in
                     ["gene_id", "transcript_id", "chrom", "strand", "spans"])


GTF_LINES = [
    # g1 '+'
    "chr1\tsyn\texon\t11\t30\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\"; gene_biotype \"protein_coding\";",
    "chr1\tsyn\texon\t51\t70\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"2\"; gene_biotype \"protein_coding\";",
    "chr1\tsyn\tstart_codon\t20\t22\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\";",
    "chr1\tsyn\tstop_codon\t62\t64\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"2\";",
    # g2 '-' (5' end = rightmost)
    "chr1\tsyn\texon\t121\t150\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; exon_number \"1\"; gene_biotype \"protein_coding\";",
    "chr1\tsyn\tstart_codon\t137\t139\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; exon_number \"1\";",
    "chr1\tsyn\tstop_codon\t128\t130\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; exon_number \"1\";",
    # g3 '+' noncoding
    "chr1\tsyn\texon\t201\t215\t.\t+\t.\tgene_id \"g3\"; transcript_id \"t3\"; exon_number \"1\"; gene_biotype \"lincRNA\";",
]


def write_synthetic(data_dir):
    """Write CHR1.fa and CHR1.gtf into ``data_dir``; return their paths."""
    import os

    fa = os.path.join(data_dir, "CHR1.fa")
    with open(fa, "w") as f:
        f.write(f">chr1\n{CHR1}\n")
    # faidx so pysam.FastaFile can fetch
    import pysam

    pysam.faidx(fa)
    gt = os.path.join(data_dir, "CHR1.gtf")
    with open(gt, "w") as f:
        f.write("\n".join(GTF_LINES) + "\n")
    return fa, gt
