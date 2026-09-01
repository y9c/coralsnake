"""liftover --table: convert a sites table between transcript(gene) and
genome coordinates (strand-aware, exon-spliced), reusing the `prepare`
annotation map (gene -> Transcript) from :mod:`coralsnake.utils`.

The mapping is the exact inverse of each other:

  t2g  (gene, genePos 1-based)        -> genome (Chrom, Pos 1-based)
  g2t  genome (Chrom, Pos, Strand)    -> (gene, genePos 1-based)

Implementation is vectorized over the flat exon table (numpy).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from .utils import load_annotation

logger = logging.getLogger(__name__)


def _flat_exons(annotation, pick_transcript="first"):
    """Flatten gene->transcripts->exons into numpy exon table.

    Returns (exons, per_gene) where `exons` is a flat rows array
    [gene, chrom, strand(+/-), tx_start, tx_end, g_start, g_end]
    and `per_gene` maps gene -> row slice (start, stop).
    """
    rows = []
    genes = []
    gene_order = []
    for gene, txmap in annotation.items():
        tx = (sorted(txmap.values(),
                     key=lambda t: (t.gene_id, t.transcript_id))
              if pick_transcript == "first" else next(iter(txmap.values())))
        if not isinstance(tx, list):
            tx = [tx]
        txs = tx if isinstance(tx, list) else [tx]
        for t in txs:
            strand = 1 if t.strand == "+" else -1
            cum = t.cum_exon_lens
            prev = 0
            for i, exon in enumerate(t.exons_forwards):
                g_start, g_end = exon.start, exon.end
                rows.append((gene, t.chrom, strand, prev,
                             cum[i], g_start, g_end))
                prev = cum[i]
                if gene not in gene_order:
                    gene_order.append(gene)
    arr = np.array(rows, dtype=object) if rows else np.empty((0, 7), dtype=object)
    return arr, gene_order


def _group_slices(arr, gene_order):
    order = {g: i for i, g in enumerate(gene_order)}
    per_gene = {}
    starts = {}

    # group by gene (vectorized via argsort on the gene column)
    idx = np.argsort(arr[:, 0], kind="stable")
    sorted_genes = arr[idx, 0]
    uniq, uniq_start = np.unique(sorted_genes, return_index=True)
    uniq_end = np.append(uniq_start[1:], len(arr))
    for u, s, e in zip(uniq, uniq_start, uniq_end):
        per_gene[u] = idx[s:e]
    return per_gene


def transcript_to_genome(annotation, genes, tx_pos1):
    """genes: np array of gene ids; tx_pos1: 1-based transcript positions.

    Returns (chroms, genome_pos1): str/None, int/0 for rows without a matching exon.
    """
    genes = np.asarray(genes)
    tx_pos1 = np.asarray(tx_pos1, dtype=np.int64)
    chroms = np.full(len(genes), None, dtype=object)
    gpos = np.zeros(len(genes), dtype=np.int64)
    for gene in np.unique(genes):
        sel = np.flatnonzero(genes == gene)
        tx = annotation.get(gene)
        if not tx:
            continue
        t = next(iter(tx.values()))
        cum = np.array(t.cum_exon_lens, dtype=np.int64)
        starts = cum - np.array([e.end - e.start for e in t.exons_forwards],
                                dtype=np.int64)
        p0 = tx_pos1[sel] - 1
        exon_i = np.searchsorted(cum, p0, side="right")  # first cum > p0
        ok = (exon_i > 0) & (p0 >= np.minimum.accumulate(starts[:-1] if len(starts) > 1 else starts))
        # simpler: exon where starts[i] <= p0 < cum[i]
        exon_i = np.searchsorted(cum, p0, side="right")
        det = (exon_i >= 1) & (exon_i <= len(cum))
        det &= (p0 >= np.take(starts, np.clip(exon_i - 1, 0, len(starts) - 1))) \
               & (p0 < np.take(cum, np.clip(exon_i - 1, 0, len(cum) - 1)))
        for k in np.flatnonzero(det):
            i = int(exon_i[k]) - 1
            exon = t.exons_forwards[i]
            offset = int(p0[k]) - int(starts[i])
            if t.strand == "+":
                gp = exon.start + offset + 1
            else:
                gp = exon.end - offset
            chroms[sel[k]] = t.chrom
            gpos[sel[k]] = gp
    return chroms, gpos


def genome_to_transcript(annotation, chroms, genome_pos1, strands):
    """genome -> (genes, tx_pos1); first matching gene wins."""
    chroms = np.asarray(chroms)
    genome_pos1 = np.asarray(genome_pos1, dtype=np.int64)
    strands = np.asarray(strands)
    genes = np.full(len(chroms), None, dtype=object)
    gpos = np.zeros(len(chroms), dtype=np.int64)

    # index exons by chrom
    chrom_index = {}
    for gene, txmap in annotation.items():
        for t in txmap.values():
            for i, exon in enumerate(t.exons_forwards):
                chrom_index.setdefault(t.chrom, []).append(
                    (gene, t, i, exon.start, exon.end))
    for ch in np.unique(chroms):
        if ch not in chrom_index:
            continue
        rows = chrom_index[ch]
        gs = np.array([r[3] for r in rows], dtype=np.int64)
        ge = np.array([r[4] for r in rows], dtype=np.int64)
        sel = np.flatnonzero(chroms == ch)
        p0 = genome_pos1[sel] - 1
        # any exon with gs <= p0 < ge
        for row in rows:
            gene, t, i, g_start, g_end = row
            det = (p0 >= g_start) & (p0 < g_end)
            if not np.any(det):
                continue
            cum = np.array(t.cum_exon_lens, dtype=np.int64)
            starts = cum - np.array([e.end - e.start for e in t.exons_forwards],
                                    dtype=np.int64)
            for k in np.flatnonzero(det):
                if genes[sel[k]] is not None:
                    continue  # first match wins
                offset = int(p0[k]) - g_start
                if t.strand == "+":
                    tp = int(starts[i]) + offset
                else:
                    tp = int(starts[i]) + (g_end - 1 - int(p0[k]))
                genes[sel[k]] = gene
                gpos[sel[k]] = tp + 1
    return genes, gpos


def run_liftover_table(input_table, output_table, annotation_file,
                       direction, gene_col="Chrom", pos_col="Pos",
                       strand_col="Strand", separator="\t"):
    """Read a table (with header), convert the position columns, write back."""
    import xopen

    annot = load_annotation(annotation_file, with_header=True)
    with xopen.xopen(input_table, "r") as f:
        header = f.readline().rstrip("\n").split(separator)
        rows = [ln.rstrip("\n").split(separator) for ln in f]
    idx = {c: i for i, c in enumerate(header)}
    if pos_col not in idx:
        raise ValueError("input table missing column %r" % pos_col)

    pos_arr = np.array([float(r[idx[pos_col]]) for r in rows], dtype=np.int64)
    out_cols = list(header)

    if direction == "t2g":
        if gene_col not in idx:
            raise ValueError("t2g input needs a %r (gene) column" % gene_col)
        genes = np.array([r[idx[gene_col]] for r in rows], dtype=object)
        chroms, gpos = transcript_to_genome(annot, genes, pos_arr)
        out_cols += ["GenomeChrom", "GenomePos"]
    elif direction == "g2t":
        chroms = np.array([r[idx[gene_col]] for r in rows], dtype=object)
        strands = np.array([r[idx[strand_col]] for r in rows], dtype=object)
        genes, gpos = genome_to_transcript(annot, chroms, pos_arr, strands)
        out_cols += ["Gene", "GenePos"]
    else:
        raise ValueError("direction must be t2g or g2t")

    with open(output_table, "w") as f:
        f.write(separator.join(out_cols) + "\n")
        for i, r in enumerate(rows):
            line = list(r)
            if direction == "t2g":
                line += [chroms[i] if chroms[i] is not None else "",
                         str(gpos[i]) if chroms[i] is not None else ""]
            else:
                line += [genes[i] if genes[i] is not None else "",
                         str(gpos[i]) if genes[i] is not None else ""]
            f.write(separator.join(str(x) for x in line) + "\n")
    logger.info("liftover table %s -> %s (%d rows)", direction, output_table, len(rows))
