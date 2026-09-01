"""liftover --table: convert a sites table between transcript(gene) and genome
coordinates (strand-aware, exon-spliced), reusing the `prepare` annotation map
(gene -> Transcript) from :mod:`coralsnake.utils`.

  t2g  (gene, genePos 1-based)     -> (GenomeChrom, GenomePos 1-based)
  g2t  (GenomeChrom, Pos, Strand)  -> (Gene, GenePos 1-based)

Both directions use transcript-order exons (5'->3'): for a '+' gene that is the
genomic-sorted exons; for a '-' gene it is their reverse.  Vectorized per gene.
"""

from __future__ import annotations

import logging

import numpy as np

from .utils import load_annotation

logger = logging.getLogger(__name__)


def _tx_exons(t):
    """(exons_in_tx_order, strand_factor)."""
    order = list(t.exons_forwards)
    if t.strand == "-":
        order.reverse()
    return order


def transcript_to_genome(annotation, genes, tx_pos1):
    """(gene, 1-based tx pos) -> (genome chroms, 1-based genome pos)."""
    genes = np.asarray(genes)
    tx_pos1 = np.asarray(tx_pos1, dtype=np.int64)
    chroms = np.full(len(genes), None, dtype=object)
    gpos = np.zeros(len(genes), dtype=np.int64)

    for gene in np.unique(genes):
        tx = annotation.get(gene)
        if not tx:
            continue
        t = next(iter(tx.values()))
        exons = _tx_exons(t)
        starts, ends = [], []
        s = 0
        for e in exons:
            starts.append(s)
            s += e.end - e.start
            ends.append(s)
        starts = np.array(starts, dtype=np.int64)
        ends = np.array(ends, dtype=np.int64)
        sel = np.flatnonzero(genes == gene)
        p0 = tx_pos1[sel] - 1
        i = np.searchsorted(ends, p0, side="right")
        ok = (i >= 0) & (i < len(starts)) & (p0 >= starts[i]) & (p0 < ends[i])
        for k in np.flatnonzero(ok):
            e = exons[int(i[k])]
            off = int(p0[k]) - int(starts[int(i[k])])
            chroms[sel[k]] = t.chrom
            gpos[sel[k]] = (e.start + off + 1) if t.strand == "+" else (e.end - off)
    return chroms, gpos


def genome_to_transcript(annotation, chroms, genome_pos1, strands):
    """(chrome, 1-based genome pos, strand) -> (gene, 1-based tx pos).

    First exon (in tx order) that contains the position wins per site.
    """
    chroms = np.asarray(chroms)
    genome_pos1 = np.asarray(genome_pos1, dtype=np.int64)
    strands = np.asarray(strands)
    genes = np.full(len(chroms), None, dtype=object)
    gpos = np.zeros(len(chroms), dtype=np.int64)

    chrom_index = {}
    for gene, txmap in annotation.items():
        for t in txmap.values():
            for exon in _tx_exons(t):
                chrom_index.setdefault(t.chrom, []).append((gene, t, exon))

    for ch in set(c for c in chroms.tolist() if c is not None):
        cand = chrom_index.get(ch)
        if not cand:
            continue
        sel = np.flatnonzero(chroms == ch)
        p0 = genome_pos1[sel] - 1
        for gene, t, exon in cand:
            det = (p0 >= exon.start) & (p0 < exon.end)
            if not np.any(det):
                continue
            exons = _tx_exons(t)
            starts = []
            s = 0
            for e in exons:
                starts.append(s)
                s += e.end - e.start
            for k in np.flatnonzero(det):
                if genes[sel[k]] is not None:
                    continue  # first match wins
                # position within the transcript (which exon index holds p0)
                off = int(p0[k]) - exon.start
                tp = starts[exons.index(exon)] + (off if t.strand == "+" else
                                                   (exon.end - 1 - int(p0[k])))
                genes[sel[k]] = gene
                gpos[sel[k]] = tp + 1
        # keep per-site detection loop; flag not found -> empty
    return genes, gpos


def run_liftover_table(input_table, output_table, annotation_file,
                       direction, gene_col="Chrom", pos_col="Pos",
                       strand_col="Strand", separator="\t"):
    """Read a table (with header), convert positions, write back."""
    import xopen

    with xopen.xopen(input_table, "r") as f:
        header = f.readline().rstrip("\n").split(separator)
        rows = [ln.rstrip("\n").split(separator) for ln in f]
    idx = {c: i for i, c in enumerate(header)}
    if pos_col not in idx:
        raise ValueError("input table missing column %r" % pos_col)
    pos_arr = np.array([float(r[idx[pos_col]]) for r in rows], dtype=np.int64)

    annot = load_annotation(annotation_file, with_header=True)
    out_cols = list(header)

    if direction == "t2g":
        if gene_col not in idx:
            raise ValueError("t2g input needs a %r (gene) column" % gene_col)
        genes = np.array([r[idx[gene_col]] for r in rows], dtype=object)
        chroms, gpos = transcript_to_genome(annot, genes, pos_arr)
        out_cols += ["GenomeChrom", "GenomePos"]
    elif direction == "g2t":
        if gene_col not in idx or strand_col not in idx:
            raise ValueError("g2t input needs %r and %r columns" % (gene_col, strand_col))
        chr_arr = np.array([r[idx[gene_col]] for r in rows], dtype=object)
        st_arr = np.array([r[idx[strand_col]] for r in rows], dtype=object)
        genes, gpos = genome_to_transcript(annot, chr_arr, pos_arr, st_arr)
        out_cols += ["Gene", "GenePos"]
    else:
        raise ValueError("direction must be t2g or g2t")

    with open(output_table, "w") as f:
        f.write(separator.join(out_cols) + "\n")
        for i, r in enumerate(rows):
            line = list(r)
            hit = (chroms[i] if direction == "t2g" else genes[i]) is not None
            if direction == "t2g":
                line += [chroms[i] if hit else "", str(gpos[i]) if hit else ""]
            else:
                line += [genes[i] if hit else "", str(gpos[i]) if hit else ""]
            f.write(separator.join(str(x) for x in line) + "\n")
    logger.info("liftover table %s -> %s (%d rows)", direction, output_table, len(rows))
