#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# RNA-seq quality control for coralsnake.
#
# This is deliberately a *slight* extension of coralsnake: rather than
# re-implementing annotation / interval machinery, it leans on what the package
# already ships:
#
#   * coralsnake.genemodel.GeneModel      -> parse the GTF into genes + exons
#   * coralsnake.utils.interval_groups + ruranges.numpy.overlaps
#                                         -> vectorized per-read interval overlap
#   * pysam.AlignmentFile                 -> stream the BAM
#   * xopen                               -> gzip / stdin-stdout aware text output
#
# It computes the core QC metrics that need only a genome-aligned BAM + a GTF
# (no optional BED/FASTA): the mapping-quality statistics, the exonic / intronic
# / intergenic / ambiguous / intragenic read classification, the rRNA rate, gene
# read counts, TPM, genes detected, 3'/5' bias and fragment-size statistics (the
# latter two are derived from the GTF exons + reads, so no separate BED is
# needed).  Per-gene/per-exon base coverage + CV and GC content (which need a
# FASTA) and the legacy counting rules are left out.
#
# Metric definitions follow the standard RNA-seq QC metrics (see the output
# documentation; the classification rule is the same used by popular RNA-seq
# QC tools: a read is exonic when its exons fully cover every CIGAR block of a
# single gene, ambiguous when exons are hit but no single gene covers all
# blocks, intronic when it hits a gene body but no exon, intergenic otherwise).
#
# Coordinates are 0-based half-open internally (the coralsnake convention),
# converted to/from GTF 1-based form at the boundary in :func:`_build_features`.

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam
from ruranges.numpy import overlaps

from .genemodel import GeneModel
from .utils import get_logger, interval_groups

LOGGER = get_logger(__name__)

# ---------------------------------------------------------------------------
# BAM flag bits
# ---------------------------------------------------------------------------
SECONDARY = 0x100
SUPPLEMENTARY = 0x800
QCFAIL = 0x200
DUPLICATE = 0x400
PAIRED = 0x1
FIRST_IN_PAIR = 0x40
REVERSE = 0x10
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
MATE_REVERSE = 0x20
PROPER_PAIR = 0x2
SECOND_IN_PAIR = 0x80

# Genes treated as "globin" for the non-globin duplicate rate.
BLACKLISTED_GLOBINS = frozenset(
    {
        "HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBM",
        "HBQ1", "HBZ", "HBBP1", "HBZP1",
    }
)


def _is_ribosomal(*texts) -> bool:
    """Heuristic for an annotated rRNA gene.

    The standard rule matches ``transcript_type`` against the regex ``rRNA``.  Real GTFs
    vary (``gene_type``/``biotype``/``gene_name``), so we accept any of those
    attributes containing ``rrna`` / ``ribosomal`` (case-insensitive).
    """
    for t in texts:
        if not t:
            continue
        tl = str(t).lower()
        if "rrna" in tl or "ribosomal" in tl:
            return True
    return False


def _block_label(chrom: str, strand: str, stranded: bool) -> str:
    """Group label for ruranges: chromosome (plus strand when ``stranded``)."""
    if stranded:
        return chrom + ("+" if strand == "+" else "-")
    return chrom


# ---------------------------------------------------------------------------
# Feature (GTF) extraction
# ---------------------------------------------------------------------------


def _build_features(gtf: str, stranded: bool):
    """Build exon- and gene-body interval arrays from a GTF via ``GeneModel``.

    Returns ``(arrays, gene_info)`` where ``arrays`` holds parallel arrays for
    the exon and gene intervals and ``gene_info`` maps ``gene_id -> (name, length,
    strand, seqname, ribosomal)`` (length = summed exon length, used for TPM).

    Intervals are 0-based half-open; ``label`` is the ruranges group label.
    """
    model = GeneModel(gtf)

    e_start, e_end, e_label, e_strand, e_gene, e_name, e_xid, e_ribo = (
        [], [], [], [], [], [], [], []
    )
    gene_span: dict[str, tuple[int, int]] = {}   # gene row, 0-based half-open
    gene_ribo: dict[str, bool] = {}
    gene_len: dict[str, int] = defaultdict(int)
    gene_name: dict[str, str] = {}
    exon_bbox: dict[str, list[int]] = defaultdict(lambda: [None, None])
    gene_strand: dict[str, str] = {}
    gene_seqname: dict[str, str] = {}

    for row in model.iter_rows():
        gid = row.attributes.get("gene_id")
        if gid is None:
            continue
        r_name = row.attributes.get("gene_name", gid)
        gene_name.setdefault(gid, r_name)
        ribo = _is_ribosomal(
            row.attributes.get("transcript_type"),
            row.attributes.get("gene_type"),
            row.attributes.get("biotype"),
            row.attributes.get("gene_name"),
        )
        if ribo:
            gene_ribo[gid] = True
        gene_strand.setdefault(gid, row.strand)
        gene_seqname.setdefault(gid, row.seqname)

        if row.feature == "gene":
            gene_span[gid] = (row.start - 1, row.end)
        elif row.feature == "exon":
            s0, e0 = row.span_0
            b = exon_bbox[gid]
            b[0] = s0 if b[0] is None else min(b[0], s0)
            b[1] = e0 if b[1] is None else max(b[1], e0)
            gene_len[gid] += e0 - s0
            e_start.append(s0)
            e_end.append(e0)
            e_strand.append(row.strand)
            e_label.append(_block_label(row.seqname, row.strand, stranded))
            e_gene.append(gid)
            e_name.append(r_name)
            e_xid.append(row.attributes.get("exon_id") or f"{gid}_{len(e_xid)}")
            e_ribo.append(ribo)

    # gene-body intervals: gene row if present, else the exon bounding box
    g_start, g_end, g_label, g_strand, g_gene, g_name, g_ribo = (
        [], [], [], [], [], [], []
    )
    gene_info = {}
    for gid, (s0, e0) in exon_bbox.items():
        if s0 is None:
            continue
        if gid in gene_span:
            s0, e0 = gene_span[gid]
        g_start.append(s0)
        g_end.append(e0)
        strand = gene_strand.get(gid, ".")
        seqname = gene_seqname.get(gid, "")
        g_strand.append(strand)
        g_label.append(_block_label(seqname, strand, stranded))
        g_gene.append(gid)
        g_name.append(gene_name.get(gid, gid))
        g_ribo.append(gene_ribo.get(gid, False))
        gene_info[gid] = (g_name[-1], int(gene_len.get(gid, 0)), strand,
                          seqname, gene_ribo.get(gid, False))

    def asarr(x, dtype):
        return np.asarray(x, dtype=dtype) if len(x) else np.empty(0, dtype=dtype)

    # per-gene exon spans (genomic order) + each exon's rank within its gene
    gene_idx_lists = defaultdict(list)
    for i, gid in enumerate(e_gene):
        gene_idx_lists[gid].append(i)
    exon_rank = np.zeros(len(e_gene), dtype=np.int64)
    exon_tx = np.zeros(len(e_gene), dtype=np.int64)  # 5'->3' transcript offset
    gene_exons: dict[str, list[tuple[int, int]]] = {}
    for gid, idxs in gene_idx_lists.items():
        idxs = sorted(idxs, key=lambda i: (e_start[i], e_end[i]))
        for rank, i in enumerate(idxs):
            exon_rank[i] = rank
        gene_exons[gid] = [(int(e_start[i]), int(e_end[i])) for i in idxs]
        if gene_strand.get(gid) == "-":
            running = 0
            for i in reversed(idxs):
                exon_tx[i] = running
                running += e_end[i] - e_start[i]
        else:
            running = 0
            for i in idxs:
                exon_tx[i] = running
                running += e_end[i] - e_start[i]

    arrays = {
        "exon": {
            "start": asarr(e_start, np.int64),
            "end": asarr(e_end, np.int64),
            "label": asarr(e_label, object),
            "strand": asarr(e_strand, object),
            "gene": np.asarray(e_gene, dtype=object),
            "gene_name": np.asarray(e_name, dtype=object),
            "exon_id": np.asarray(e_xid, dtype=object),
            "ribosomal": asarr(e_ribo, bool),
            "rank": exon_rank,
            "tx_offset": exon_tx,
        },
        "gene": {
            "start": asarr(g_start, np.int64),
            "end": asarr(g_end, np.int64),
            "label": asarr(g_label, object),
            "strand": asarr(g_strand, object),
            "gene": np.asarray(g_gene, dtype=object),
            "gene_name": np.asarray(g_name, dtype=object),
            "ribosomal": asarr(g_ribo, bool),
        },
    }
    return arrays, gene_info, gene_exons


def _cigar_blocks(read: pysam.AlignedSegment):
    """Reference-aligned intervals (0-based half-open) from the CIGAR.

    M/=/X advance the reference and emit a block; N/D advance the reference;
    H/P/S/I consume no reference.
    """
    cig = read.cigartuples
    if not cig:
        return []
    pos = read.reference_start
    out = []
    for op, length in cig:
        if op in (0, 7, 8):  # MATCH, EQUAL, DIFF
            out.append((pos, pos + length))
            pos += length
        elif op in (3, 2):  # N, D
            pos += length
    return out


def _rate(num, den):
    if den in (0, None):
        return float("nan")
    return num / den


def _stats(values):
    """(mean, median, std, MAD, 25th, 75th) over a list of numbers.

    std is the population std (divides by n); MAD is scaled by 1.4826.
    """
    if not values:
        return float("nan"), float("nan"), float("nan"), float("nan"), \
            float("nan"), float("nan")
    a = np.asarray(values, dtype=float)
    avg = float(a.mean())
    med = float(np.median(a))
    std = float(a.std(ddof=0))
    mad = float(np.median(np.abs(a - med))) * 1.4826
    p25 = float(np.percentile(a, 25))
    p75 = float(np.percentile(a, 75))
    return avg, med, std, mad, p25, p75


def _fragment_sizes(n_reads, b_read, fc_blk, fc_ex, ex_f, rec_paired, rec_hq,
                    rec_qname, rec_tlen, per_read_total, limit):
    """Fragment sizes (insert sizes) of high-quality read pairs whose two mates
    align fully within a single common exon.

    Some tools use a ``--bed`` of non-overlapping exons here; coralsnake uses the
    GTF exons serve the same purpose.  Returns a sorted list of sizes.
    """
    if not n_reads or not len(fc_blk) or limit <= 0:
        return []
    import polars as pl

    fc_read = b_read[fc_blk]
    df = pl.DataFrame({
        "read": fc_read,
        "exon": ex_f["exon_id"][fc_ex],
        "block": fc_blk,
    })
    grp = df.group_by("read").agg(
        n_exon=pl.col("exon").n_unique(),
        n_fc_block=pl.col("block").n_unique(),
    )
    tot = pl.DataFrame({"read": np.arange(n_reads), "total": per_read_total})
    single = (
        grp.join(tot, on="read")
        .filter((pl.col("n_exon") == 1) & (pl.col("n_fc_block") == pl.col("total")))
        .join(df, on="read")
        .unique(subset=["read"])
    )
    ok = pl.DataFrame({
        "read": np.arange(n_reads),
        "paired": rec_paired, "hq": rec_hq,
        "qname": rec_qname, "tlen": rec_tlen,
    })
    cand = (
        single.join(ok, on="read")
        .filter(pl.col("paired") & pl.col("hq"))
        .select(["qname", "exon", "tlen"])
    )
    # a fragment is sampled only when both mates of the pair share an exon
    frag = (
        cand.group_by(["qname", "exon"])
        .agg(n=pl.len(), tlen=pl.col("tlen").first())
        .filter(pl.col("n") == 2)
    )
    sizes = sorted(int(v) for v in frag["tlen"].to_list())
    return sizes[:limit]


def _three_prime_bias(b_start, b_end, b_read, fc_ex, fc_blk, ex_f, rec_hq,
                      do_exon, gene_exons, gene_info, unique_gene_count,
                      offset, window, gene_min_len, det_threshold):
    """Per-gene 3' bias = cov3 / (cov3 + cov5), from GTF exons + read coverage.

    For each eligible gene (exon-total length >= gene_min_len and >= det_threshold
    unique reads), the spliced-transcript per-base coverage is accumulated from
    the fully-covered exon blocks of high-quality unambiguous reads, then the
    median coverage of a ``window``-bp 5' and 3' window (each offset by
    ``offset`` bp into the gene) is compared.  Returns a list of biases.
    """
    eligible = {
        g for g in gene_exons
        if gene_info.get(g, (0, 0))[1] >= gene_min_len
        and unique_gene_count.get(g, 0) >= det_threshold
    }
    if not eligible or not len(fc_blk):
        return []

    # Vectorized transcript-position mapping for every fully-covered exon block.
    # keep = high-quality + unambiguous + gene eligible.
    fc_read = b_read[fc_blk]
    ex_idx = fc_ex
    keep = rec_hq[fc_read] & do_exon[fc_read]
    gene_ids = ex_f["gene"][ex_idx]
    gene_code = {g: i for i, g in enumerate(eligible)}
    code = np.fromiter(
        (gene_code.get(g, -1) for g in gene_ids), dtype=np.int64, count=len(gene_ids)
    )
    keep &= code >= 0

    pos = ex_f["tx_offset"][ex_idx] + np.where(
        ex_f["strand"][ex_idx] == "+",
        b_start[fc_blk] - ex_f["start"][ex_idx],
        ex_f["end"][ex_idx] - b_end[fc_blk],
    )
    blen = b_end[fc_blk] - b_start[fc_blk]
    codes = code[keep]
    pos = pos[keep]
    blen = blen[keep]

    # per-gene per-transcript-position coverage via a difference array:
    #   diff[pos] += 1 ; diff[pos+len] -= 1 ; coverage = cumsum(diff)
    elig_list = list(eligible)
    gene_len = [int(gene_info[g][1]) for g in elig_list]
    biases = []
    for gc in np.unique(codes):
        m = codes == gc
        p = pos[m]
        bl = blen[m]
        length = gene_len[gc]
        diff = np.zeros(length + 1, dtype=np.int64)
        np.add.at(diff, p, 1)
        np.add.at(diff, np.minimum(p + bl, length), -1)
        cov = np.cumsum(diff)[:length]
        if length < window or length <= 2 * offset:
            continue
        cov5 = float(np.median(cov[offset: offset + window]))
        cov3 = float(np.median(cov[length - offset - window: length - offset]))
        if cov5 + cov3 <= 0:
            continue
        biases.append(cov3 / (cov5 + cov3))
    return biases


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def run_rnaseqc(
    bam: str,
    gtf: str,
    outdir: str,
    sample: str | None = None,
    unpaired: bool = False,
    mapping_quality: int = 255,
    base_mismatch: int = 6,
    detection_threshold: int = 5,
    stranded: str | None = None,
    write_counts: bool = True,
    bias_offset: int = 150,
    bias_window: int = 100,
    bias_gene_length: int = 600,
    fragment_samples: int = 1_000_000,
):
    """Run the core RNA-seq QC metrics on a genome-aligned BAM + GTF.

    Parameters
    ----------
    bam : input (coordinate-sorted) SAM/BAM/CRAM.
    gtf : gene annotation (GTF/GFF); contig names must match the BAM.
    outdir : directory to write ``{sample}.metrics.tsv`` (and, when
        ``write_counts``, ``{sample}.gene_reads.tsv`` / ``.gene_tpm.tsv`` /
        ``.exon_reads.tsv``).
    sample : sample name; default: the BAM filename (without extension).
    unpaired : allow single-end libraries (do not require proper pairs).
    mapping_quality : minimum mapq for a read to be "high quality".
    base_mismatch : max NM mismatches for a read to be "high quality".
    detection_threshold : gene reads (unique) required for "gene detected".
    stranded : optional ``"RF"/"FR"/"rf"/"fr"``; restricts features to the
        read's strand (strand-specific libraries).
    write_counts : also emit the gene/exon count tables.
    bias_offset : 3' bias offset into the gene (bp) from each end.
    bias_window : 3' bias window size (bp) at each end.
    bias_gene_length : minimum exon-total gene length (bp) for 3' bias.
    fragment_samples : max number of fragment-size samples to collect (from
        pairs aligning to the same exon).

    Returns
    -------
    dict of ``{metric_name: value}`` (also written to ``metrics.tsv``).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if sample is None:
        sample = Path(bam).name
        for suf in (".bam", ".sam", ".cram"):
            if sample.lower().endswith(suf):
                sample = sample[: -len(suf)]
                break

    if stranded:
        s = str(stranded).lower()
        if s not in ("rf", "fr"):
            raise ValueError("--stranded must be in {'RF','RF','rf','fr'}")
    use_strand = stranded is not None

    arrays, gene_info, gene_exons = _build_features(gtf, stranded=use_strand)
    ex_f, g_f = arrays["exon"], arrays["gene"]

    mode = "rb" if bam.lower().endswith((".bam", ".cram")) else "r"
    bam_handle = pysam.AlignmentFile(bam, mode)

    # BAM contigs must overlap the annotation.
    bam_refs = set(bam_handle.references)
    gtf_refs = set(ex_f["label"]) | set(g_f["label"])
    if use_strand:
        gtf_refs = {lab[:-1] for lab in gtf_refs}
    else:
        gtf_refs = set(lab for lab in ex_f["label"]) | set(lab for lab in g_f["label"])
    if not (bam_refs & gtf_refs):
        bam_handle.close()
        raise ValueError("BAM shares no contigs with the GTF annotation.")

    # ---- mapping-level counters -------------------------------------------
    c = {n: 0 for n in (
        "total", "secondary", "supplementary", "qcfail", "low_mapq",
        "unique_pass", "unpaired_reads", "mapped", "mapped_dup",
        "mapped_unique", "end1_mapped", "end2_mapped", "end1_mism",
        "end2_mism", "end1_bases", "end2_bases", "mismatched_bases",
        "total_bases", "total_pairs", "dup_pairs", "unique_frags",
        "high_q", "low_q", "reads_used", "alignment_blocks",
    )}

    # per-read records for the deferred annotation classification
    rec_aligned: list[int] = []   # aligned bases (M/=,X)
    rec_hq: list[bool] = []
    rec_dup: list[bool] = []
    rec_first: list[bool] = []
    rec_reverse: list[bool] = []
    rec_paired: list[bool] = []
    rec_qname: list[str] = []      # for fragment-size pairing
    rec_tlen: list[int] = []       # abs(insert size / TLEN)
    rec_pos_end: list[int] = []    # alignment end (0-based half-open)
    rec_mate_reverse: list[bool] = []
    rec_mate_pos: list[int] = []
    b_start: list[int] = []
    b_end: list[int] = []
    b_read: list[int] = []
    b_label: list[str] = []

    read_length = 0

    for read in bam_handle:
        flag = read.flag
        c["total"] += 1

        if flag & SECONDARY:
            c["secondary"] += 1
        if flag & SUPPLEMENTARY:
            c["supplementary"] += 1
        elif flag & QCFAIL:
            c["qcfail"] += 1
        elif read.mapping_quality < mapping_quality:
            c["low_mapq"] += 1

        if flag & (SECONDARY | QCFAIL | SUPPLEMENTARY):
            continue  # not in "Unique Mapping, Vendor QC Passed Reads"

        c["unique_pass"] += 1
        if not (flag & PAIRED):
            c["unpaired_reads"] += 1
        if flag & UNMAPPED:
            continue

        c["mapped"] += 1
        if flag & DUPLICATE:
            c["mapped_dup"] += 1
        else:
            c["mapped_unique"] += 1

        blocks = _cigar_blocks(read)
        aligned = sum(e - s for s, e in blocks)
        read_length = max(read_length, aligned)
        nm = read.get_tag("NM") if read.has_tag("NM") else 0

        if flag & PAIRED and not (flag & MATE_UNMAPPED):
            if flag & FIRST_IN_PAIR:
                c["total_pairs"] += 1
                c["end1_mapped"] += 1
                c["end1_mism"] += nm
                c["end1_bases"] += aligned
                c["dup_pairs"] += 1 if flag & DUPLICATE else 0
                if not (flag & DUPLICATE):
                    c["unique_frags"] += 1
            else:
                c["end2_mapped"] += 1
                c["end2_mism"] += nm
                c["end2_bases"] += aligned
        c["total_bases"] += aligned
        c["mismatched_bases"] += nm

        high_q = (
            nm <= base_mismatch
            and (unpaired or bool(flag & PROPER_PAIR))
            and read.mapping_quality >= mapping_quality
        )
        c["high_q" if high_q else "low_q"] += 1
        c["reads_used"] += 1
        c["alignment_blocks"] += len(blocks)

        if not blocks or read.reference_name not in gtf_refs:
            continue

        ridx = len(rec_aligned)
        rec_aligned.append(aligned)
        rec_hq.append(high_q)
        rec_dup.append(bool(flag & DUPLICATE))
        rec_first.append(bool(flag & FIRST_IN_PAIR) or unpaired)
        rec_reverse.append(bool(flag & REVERSE))
        rec_paired.append(bool(flag & PAIRED))
        rec_qname.append(read.query_name)
        rec_tlen.append(abs(read.template_length) if flag & PAIRED else 0)
        rec_pos_end.append(read.reference_end)
        rec_mate_reverse.append(bool(flag & MATE_REVERSE))
        rec_mate_pos.append(read.next_reference_start if flag & PAIRED else -1)
        b_read.extend([ridx] * len(blocks))
        rstrand = "-" if flag & REVERSE else "+"
        b_label.extend(
            _block_label(read.reference_name, rstrand, use_strand) for _ in blocks
        )
        b_start.extend(s for s, _ in blocks)
        b_end.extend(e for _, e in blocks)

    bam_handle.close()

    n_reads = len(rec_aligned)
    if len(b_start):
        b_start = np.asarray(b_start, dtype=np.int64)
        b_end = np.asarray(b_end, dtype=np.int64)
        b_read = np.asarray(b_read, dtype=np.int64)
        b_label = np.asarray(b_label, dtype=object)
        b_len = b_end - b_start

    rec_hq = np.asarray(rec_hq, dtype=bool)
    rec_dup = np.asarray(rec_dup, dtype=bool)
    rec_first = np.asarray(rec_first, dtype=bool)
    rec_reverse = np.asarray(rec_reverse, dtype=bool)
    rec_paired = np.asarray(rec_paired, dtype=bool)
    rec_aligned = np.asarray(rec_aligned, dtype=np.int64)
    rec_qname = np.asarray(rec_qname, dtype=object)
    rec_tlen = np.asarray(rec_tlen, dtype=np.int64)
    rec_pos_end = np.asarray(rec_pos_end, dtype=np.int64)
    rec_mate_reverse = np.asarray(rec_mate_reverse, dtype=bool)
    rec_mate_pos = np.asarray(rec_mate_pos, dtype=np.int64)

    # ---- annotation classification (vectorized) ---------------------------
    exonic_any = np.zeros(n_reads, dtype=bool)   # any exon overlap
    intragenic = np.zeros(n_reads, dtype=bool)   # any gene-body overlap
    plus = np.zeros(n_reads, dtype=bool)
    minus = np.zeros(n_reads, dtype=bool)
    ribosomal = np.zeros(n_reads, dtype=bool)

    gene_count: dict[str, int] = defaultdict(int)
    unique_gene_count: dict[str, int] = defaultdict(int)
    exon_count: dict[str, float] = defaultdict(float)
    do_exon = np.zeros(n_reads, dtype=bool)
    non_globin = non_globin_dup = 0

    # fully-covered (exon, block) pairs, lifted to function scope so the 3'
    # bias and fragment-size metrics can reuse them after classification
    fc_ex = np.empty(0, dtype=np.int64)
    fc_blk = np.empty(0, dtype=np.int64)

    per_read_total = np.bincount(b_read, minlength=n_reads)

    if n_reads and len(b_start):
        def _overlap(feats):
            eg, bg = interval_groups(feats["label"], b_label)
            if not len(eg):
                return np.empty(0, np.int64), np.empty(0, np.int64)
            return overlaps(
                starts=feats["start"], ends=feats["end"],
                starts2=b_start, ends2=b_end,
                groups=eg.astype(np.uint32), groups2=bg.astype(np.uint32),
            )

        # exon overlaps
        if len(ex_f["start"]):
            ex_idx, blk = _overlap(ex_f)
            if len(ex_idx):
                hit = b_read[blk]
                exonic_any[hit] = True
                plus[hit] |= ex_f["strand"][ex_idx] == "+"
                minus[hit] |= ex_f["strand"][ex_idx] == "-"
                ribosomal[hit] |= ex_f["ribosomal"][ex_idx]

                # a block is fully exonic when an exon covers it entirely
                fc = (ex_f["start"][ex_idx] <= b_start[blk]) & (
                    ex_f["end"][ex_idx] >= b_end[blk]
                )
                if fc.any():
                    fc_ex = ex_idx[fc]
                    fc_blk = blk[fc]
                    fc_read = b_read[fc_blk]
                    fc_gene = ex_f["gene"][fc_ex]

                    import polars as pl

                    fc_df = pl.DataFrame({
                        "read": fc_read,
                        "block": fc_blk,
                        "gene": fc_gene,
                        "exon": fc_ex,
                    })
                    # unambiguous = a gene whose exons fully cover EVERY block
                    unamb = (
                        fc_df.group_by(["read", "gene"])
                        .agg(n_fc=pl.col("block").n_unique())
                        .join(
                            pl.DataFrame({
                                "read": np.arange(n_reads),
                                "total": per_read_total,
                            }),
                            on="read",
                        )
                        .filter(pl.col("n_fc") == pl.col("total"))
                    )
                    do_exon_reads = unamb["read"].unique().to_numpy()
                    do_exon[do_exon_reads] = True

                    # gene counts --- high-quality reads only
                    hq = pl.DataFrame({
                        "read": np.arange(n_reads), "hq": rec_hq, "dup": rec_dup,
                    })
                    unamb_hq = unamb.join(hq, on="read").filter(pl.col("hq"))
                    for g, cnt in unamb_hq.group_by("gene").len().iter_rows():
                        gene_count[g] += cnt
                    for g, cnt in (
                        unamb_hq.filter(~pl.col("dup")).group_by("gene").len()
                        .iter_rows()
                    ):
                        unique_gene_count[g] += cnt

                    # non-globin --- vectorized, high-quality unambiguous reads
                    hq_unamb = unamb_hq.select("read")
                    hq_unamb_set = set(hq_unamb["read"].to_list())
                    if hq_unamb.height:
                        gn = np.array(
                            [gene_info.get(g, (g,))[0]
                             for g in unamb_hq["gene"].to_list()],
                            dtype=object,
                        )
                        per_read_globin = (
                            pl.DataFrame({
                                "read": unamb_hq["read"],
                                "globin": np.array(
                                    [x in BLACKLISTED_GLOBINS for x in gn]
                                ),
                            })
                            .group_by("read")
                            .agg(any_globin=pl.col("globin").any())
                        )
                        nong = per_read_globin.filter(~pl.col("any_globin"))
                        non_globin = nong.height
                        non_globin_dup = int(
                            nong.join(
                                pl.DataFrame({"read": np.arange(n_reads),
                                              "dup": rec_dup}),
                                on="read",
                            ).filter(pl.col("dup")).height
                        )

                    dose = fc_df.filter(
                        pl.col("read").is_in(hq_unamb_set)
                    )
                    blk_arr = dose["block"].to_numpy()
                    dose = dose.with_columns(
                        alen=pl.Series(rec_aligned[dose["read"].to_numpy()])
                    ).with_columns(
                        d=pl.Series(b_len[blk_arr]) / pl.col("alen")
                    )
                    for ex, d in dose.group_by("exon").agg(
                        pl.col("d").sum()
                    ).iter_rows():
                        exon_count[ex_f["exon_id"][ex]] += d

        # gene-body overlaps (intragenic) + rRNA
        if len(g_f["start"]):
            gi, blk2 = _overlap(g_f)
            if len(gi):
                hit2 = b_read[blk2]
                intragenic[hit2] = True
                plus[hit2] |= g_f["strand"][gi] == "+"
                minus[hit2] |= g_f["strand"][gi] == "-"
                ribosomal[hit2] |= g_f["ribosomal"][gi]

    # ---- classification counts -------------------------------------------
    exonic_ct = int(do_exon.sum())
    intronic_ct = int((~exonic_any & intragenic).sum())
    intergenic_ct = int((~exonic_any & ~intragenic).sum())
    ambiguous_ct = int((exonic_any & ~do_exon).sum())
    intragenic_ct = exonic_ct + intronic_ct

    end1_sense = end1_antisense = end2_sense = end2_antisense = 0
    if n_reads:
        sense_ok = (plus ^ minus) & (rec_paired | True)
        for i in np.nonzero(sense_ok)[0]:
            sense = (rec_reverse[i] and minus[i]) or (
                (not rec_reverse[i]) and plus[i]
            )
            if rec_first[i]:
                if sense:
                    end1_sense += 1
                else:
                    end1_antisense += 1
            else:
                if sense:
                    end2_sense += 1
                else:
                    end2_antisense += 1

    # ---- fragment sizes & 3' bias (from GTF exons + reads) ----------------
    frag_sizes = _fragment_sizes(
        n_reads, b_read, fc_blk, fc_ex, ex_f, rec_paired, rec_hq, rec_qname,
        rec_tlen, per_read_total, fragment_samples)
    biases = _three_prime_bias(
        b_start, b_end, b_read, fc_ex, fc_blk, ex_f, rec_hq, do_exon,
        gene_exons, gene_info, unique_gene_count, bias_offset, bias_window,
        bias_gene_length, detection_threshold)

    unique_pass = c["unique_pass"]
    mapped = c["mapped"]

    metrics = {}
    metrics["Sample"] = sample
    metrics["Mapping Rate"] = _rate(c["mapped"], unique_pass)
    metrics["Unique Rate of Mapped"] = _rate(c["mapped_unique"], mapped)
    metrics["Duplicate Rate of Mapped"] = _rate(c["mapped_dup"], mapped)
    metrics["Duplicate Rate of Mapped, excluding Globins"] = _rate(
        non_globin_dup, non_globin)
    metrics["Base Mismatch"] = _rate(c["mismatched_bases"], c["total_bases"])
    metrics["End 1 Mapping Rate"] = 2.0 * _rate(c["end1_mapped"], unique_pass)
    metrics["End 2 Mapping Rate"] = 2.0 * _rate(c["end2_mapped"], unique_pass)
    metrics["End 1 Mismatch Rate"] = _rate(c["end1_mism"], c["end1_bases"])
    metrics["End 2 Mismatch Rate"] = _rate(c["end2_mism"], c["end2_bases"])
    metrics["Expression Profiling Efficiency"] = _rate(exonic_ct, unique_pass)
    metrics["High Quality Rate"] = _rate(c["high_q"], mapped)
    metrics["Exonic Rate"] = _rate(exonic_ct, mapped)
    metrics["Intronic Rate"] = _rate(intronic_ct, mapped)
    metrics["Intergenic Rate"] = _rate(intergenic_ct, mapped)
    metrics["Intragenic Rate"] = _rate(intragenic_ct, mapped)
    metrics["Ambiguous Alignment Rate"] = _rate(ambiguous_ct, mapped)
    metrics["High Quality Exonic Rate"] = _rate(
        int((do_exon & rec_hq).sum()) if n_reads else 0, c["high_q"])
    metrics["High Quality Intronic Rate"] = _rate(
        int((~exonic_any & intragenic & rec_hq).sum()) if n_reads else 0,
        c["high_q"])
    metrics["High Quality Intergenic Rate"] = _rate(
        int((~exonic_any & ~intragenic & rec_hq).sum()) if n_reads else 0,
        c["high_q"])
    metrics["High Quality Intragenic Rate"] = _rate(
        int(((do_exon | (~exonic_any & intragenic)) & rec_hq).sum())
        if n_reads else 0, c["high_q"])
    metrics["High Quality Ambiguous Alignment Rate"] = _rate(
        int((exonic_any & ~do_exon & rec_hq).sum()) if n_reads else 0,
        c["high_q"])
    metrics["Discard Rate"] = _rate(c["mapped"] - c["reads_used"], mapped)
    metrics["rRNA Rate"] = _rate(int(ribosomal.sum()), mapped)
    metrics["End 1 Sense Rate"] = _rate(end1_sense, end1_sense + end1_antisense)
    metrics["End 2 Sense Rate"] = _rate(end2_sense, end2_sense + end2_antisense)
    metrics["Avg. Splits per Read"] = (
        _rate(c["alignment_blocks"], mapped) - 1.0 if mapped else float("nan")
    )

    metrics["Read Length"] = read_length
    detected = sum(1 for g, cnt in unique_gene_count.items()
                   if cnt >= detection_threshold)
    metrics["Genes Detected"] = detected
    metrics["Estimated Library Complexity"] = _estimate_library_complexity(
        c["unique_frags"], c["dup_pairs"])

    # 3' bias (from the GTF exons + reads, no BED needed)
    if biases:
        b_avg, b_med, b_std, b_mad, b25, b75 = _stats(biases)
        metrics["Genes used in 3' bias"] = len(biases)
        metrics["Mean 3' bias"] = b_avg
        metrics["Median 3' bias"] = b_med
        metrics["3' bias Std"] = b_std
        metrics["3' bias MAD"] = b_mad
        metrics["3' Bias, 25th Percentile"] = b25
        metrics["3' Bias, 75th Percentile"] = b75

    # fragment-size statistics (both mates align to the same GTF exon)
    if frag_sizes:
        f_avg, f_med, f_std, f_mad, _, _ = _stats(frag_sizes)
        metrics["Average Fragment Length"] = f_avg
        metrics["Fragment Length Median"] = f_med
        metrics["Fragment Length Std"] = f_std
        metrics["Fragment Length MAD"] = f_mad
        _write_fragments(outdir, sample, frag_sizes)

    metrics["Total Reads (Raw)"] = c["total"]
    metrics["Unique Mapping, Vendor QC Passed Reads (Raw)"] = unique_pass
    metrics["Mapped Reads (Raw)"] = mapped
    metrics["Mapped Unique Reads (Raw)"] = c["mapped_unique"]
    metrics["Mapped Duplicate Reads (Raw)"] = c["mapped_dup"]
    metrics["High Quality Reads (Raw)"] = c["high_q"]
    metrics["Exonic Reads (Raw)"] = exonic_ct
    metrics["Intronic Reads (Raw)"] = intronic_ct
    metrics["Intergenic Reads (Raw)"] = intergenic_ct

    _write_metrics(outdir / f"{sample}.metrics.tsv", metrics)
    LOGGER.info("RNA-seq QC metrics written to %s", outdir)

    if write_counts:
        _write_counts(outdir, sample, gene_count, unique_gene_count,
                      exon_count, gene_info)

    return metrics


# ---------------------------------------------------------------------------
# Library complexity (Landauer-Waterman)
# ---------------------------------------------------------------------------


def _estimate_library_complexity(unique: int, dup: int) -> int:
    """Estimate unique cDNA fragments (Lander-Waterman / Picard).

    Solves ``x * (1 - exp(-N/x)) = U`` (monotone in ``x``) by bisection, where
    ``N = unique + dup`` and ``U = unique``.
    """
    u, d = int(unique), int(dup)
    if u <= 0 or d <= 0:
        return 0
    n = u + d
    lo, hi = float(u), float(max(u, 1e6))
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if mid * (1.0 - np.exp(-n / mid)) > u:
            hi = mid
        else:
            lo = mid
    return int(round((lo + hi) / 2.0))


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _fmt_metric(v):
    """Render a metric value: nan/inf stay literal, integer-valued floats lose
    their decimal point, others use %.6g."""
    if isinstance(v, float):
        if v != v:  # nan
            return "nan"
        if v in (float("inf"), float("-inf")):
            return "inf" if v > 0 else "-inf"
        if v == int(v):
            return str(int(v))
        return f"{v:.6g}"
    return str(v)


def _write_metrics(path, metrics: dict):
    from xopen import xopen

    with xopen(str(path), "wt") as fh:
        for k, v in metrics.items():
            fh.write(f"{k}\t{_fmt_metric(v)}\n")


def _write_fragments(outdir, sample, sizes):
    """Write a fragment-size histogram: size\tcount."""
    from xopen import xopen

    with xopen(str(outdir / f"{sample}.fragmentSizes.txt"), "wt") as fh:
        fh.write("Fragment Size\tCount\n")
        for size, count in sorted(_counts(sizes).items()):
            fh.write(f"{size}\t{count}\n")


def _counts(values):
    out = {}
    for v in values:
        out[v] = out.get(v, 0) + 1
    return out


def _write_counts(outdir, sample, gene_count, unique_gene_count, exon_count,
                  gene_info):
    from xopen import xopen

    with xopen(str(outdir / f"{sample}.gene_reads.tsv"), "wt") as fh:
        fh.write("Name\tDescription\tCounts\n")
        for g in sorted(gene_count):
            name = gene_info.get(g, (g,))[0]
            fh.write(f"{g}\t{name}\t{int(gene_count[g])}\n")

    # TPM = (1000 * count / length), then scaled so the sum is 1e6.
    tpms = {}
    for g, cnt in gene_count.items():
        length = (gene_info.get(g, (0, 1))[1] or 1)
        tpms[g] = 1000.0 * cnt / max(1, length)
    total = sum(tpms.values()) / 1e6
    with xopen(str(outdir / f"{sample}.gene_tpm.tsv"), "wt") as fh:
        fh.write("Name\tDescription\tTPM\n")
        for g in sorted(tpms):
            name = gene_info.get(g, (g,))[0]
            fh.write(f"{g}\t{name}\t{tpms[g] / total if total else 0:.6g}\n")

    with xopen(str(outdir / f"{sample}.exon_reads.tsv"), "wt") as fh:
        fh.write("Name\tDescription\tCounts\n")
        for e in sorted(exon_count):
            fh.write(f"{e}\t\t{exon_count[e]:.6g}\n")
