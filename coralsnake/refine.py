#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Reference cleaning for coralsnake (`coralsnake refine`).
#
# Built on the shared gene-annotation model (:mod:`coralsnake.gene_annotation`):
#   * reads and writes the annotation with GeneAnnotation (attribute parsing,
#     seqname renaming/filtering, sorting, bgzip/tabix all in one place)
#   * the refined GTF keeps start/stop codon, UTR and other feature rows
#   * starts/ends of transcripts and genes are repaired (missing gene /
#     transcript / exon rows are created, overlapping exons/CDS merged)
#   * canonical transcripts are flagged (is_canonical + Ensembl_canonical tag so
#     `prepare` ranks them first); gene_biotype / transcript_biotype are written
#     so `prepare --filter-biotype / --with-biotype` work on GENCODE-style input

import gzip
import os
import re
from collections import defaultdict

import pysam

from .gene_annotation import AnnotationRow, GeneAnnotation
from .utils import get_logger

LOGGER = get_logger(__name__)


def load_seqname_mapper(path):
    """TSV: 1st column old seqname, 2nd column new seqname."""
    rename_mapper = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            row = line.strip().split("\t")
            if len(row) >= 2:
                rename_mapper[row[0]] = row[1].replace(" ", "_")
    return rename_mapper


def load_canonical_transcripts(path):
    """TSV: 1st column transcript ID (a version suffix like '.1' is ignored)."""
    canonicals = set()
    with open(path) as f:
        for line in f:
            row = line.strip().split("\t")
            if row and row[0]:
                canonicals.add(row[0])
    return canonicals


# ---------------------------------------------------------------------------
# row helpers (operate on the shared AnnotationRow)
# ---------------------------------------------------------------------------


def group_rows_by_feature(rows):
    data = defaultdict(list)
    for row in rows:
        data[row.feature].append(row)
    return data


def check_seqname_and_strand_consistency(rows):
    row1 = rows[0]
    gid = row1.attributes["gene_id"]
    for row2 in rows[1:]:
        if row1.seqname != row2.seqname:
            LOGGER.error(
                f"Gene {gid} has rows on different seqnames: "
                f"{row1.seqname} and {row2.seqname}"
            )
            return False
        if row1.strand != row2.strand:
            LOGGER.error(
                f"Gene {gid} has rows on different strands: {row1.strand} and {row2.strand}"
            )
            return False
    return True


def _bounding_row(rows, feature, drop_transcript_id=True):
    """Build a synthetic row spanning all ``rows`` (min start, max end)."""
    first = rows[0]
    attris = first.attributes.copy()
    if drop_transcript_id:
        attris.pop("transcript_id", None)
    row = AnnotationRow(
        seqname=first.seqname,
        source=first.source,
        feature=feature,
        start=min(r.start for r in rows),
        end=max(r.end for r in rows),
        score=first.score,
        strand=first.strand,
        frame=first.frame,
        attributes=attris,
    )
    return row


def create_gene_row_from_gene_rows(rows):
    gene_id = rows[0].attributes["gene_id"]
    data_features = group_rows_by_feature(rows)
    if "transcript" in data_features:
        LOGGER.info(f"Creating gene row from transcript rows for {gene_id}")
        return _bounding_row(data_features["transcript"], "gene")
    if "exon" in data_features:
        LOGGER.info(f"Creating gene row from exon rows for {gene_id}")
        return _bounding_row(data_features["exon"], "gene")
    if "CDS" in data_features:
        LOGGER.info(f"Creating gene row from CDS rows for {gene_id}")
        gene_row = _bounding_row(data_features["CDS"], "gene")
        gene_row.frame = "."
        # extend to include the stop codon
        if gene_row.strand == "-":
            gene_row.start -= 3
        else:
            gene_row.end += 3
        return gene_row
    return None


def _infer_unique(rows, label, keys, kind):
    """Collect a unique attribute value across rows; None when ambiguous."""
    vs = set()
    for row in rows:
        for k in keys:
            if k in row.attributes:
                vs.add(row.attributes[k])
    if len(vs) == 1:
        return vs.pop()
    if len(vs) > 1:
        LOGGER.info(f"{kind} {label} has multiple values for {keys}: {sorted(vs)}")
    return None


def infer_gene_type(rows):
    gid = rows[0].attributes["gene_id"]
    v = _infer_unique(rows, gid, ["gene_type", "gene_biotype"], "Gene")
    if v is not None:
        return v
    for row in rows:
        if row.feature == "CDS":
            return "protein_coding"
    return "unknown"


def infer_gene_name(rows):
    gid = rows[0].attributes["gene_id"]
    v = _infer_unique(rows, gid, ["gene_name", "gene"], "Gene")
    return v if v is not None else gid


def infer_transcript_type(rows, tid):
    v = _infer_unique(
        rows, tid, ["transcript_type", "transcript_biotype"], "Transcript"
    )
    if v is not None:
        return v
    for row in rows:
        if row.feature == "CDS":
            return "protein_coding"
    return "unknown"


def infer_transcript_name(rows, tid):
    v = _infer_unique(rows, tid, ["transcript_name"], "Transcript")
    return v if v is not None else tid


def create_transcript_row_from_transcript_rows(rows):
    data_features = group_rows_by_feature(rows)
    if "exon" in data_features:
        return _bounding_row(
            data_features["exon"], "transcript", drop_transcript_id=False
        )
    if "CDS" in data_features:
        transcript_row = _bounding_row(
            data_features["CDS"], "transcript", drop_transcript_id=False
        )
        transcript_row.frame = "."
        if transcript_row.strand == "-":
            transcript_row.start -= 3
        else:
            transcript_row.end += 3
        return transcript_row
    return None


def create_transcript_row_from_gene_row(gene_row):
    gene_attris = gene_row.attributes
    transcript_id = gene_attris["gene_id"] + "_transcript"
    transcript_attris = gene_attris.copy()
    transcript_attris["transcript_id"] = transcript_id
    transcript_attris["transcript_name"] = transcript_id
    transcript_attris["transcript_type"] = gene_attris["gene_type"]
    return AnnotationRow(
        seqname=gene_row.seqname,
        source=gene_row.source,
        feature="transcript",
        start=gene_row.start,
        end=gene_row.end,
        score=gene_row.score,
        strand=gene_row.strand,
        frame=gene_row.frame,
        attributes=transcript_attris,
    )


def create_exon_row_from_transcript_row(transcript_row):
    exon_attris = transcript_row.attributes.copy()
    exon_attris["exon_number"] = 1
    return AnnotationRow(
        seqname=transcript_row.seqname,
        source=transcript_row.source,
        feature="exon",
        start=transcript_row.start,
        end=transcript_row.end,
        score=transcript_row.score,
        strand=transcript_row.strand,
        frame=".",
        attributes=exon_attris,
    )


def create_exon_rows_from_cds_rows(cds_rows):
    exon_rows = []
    strand = cds_rows[0].strand
    for i, cds_row in enumerate(cds_rows):
        exon_attris = cds_row.attributes.copy()
        exon_row = AnnotationRow(
            seqname=cds_row.seqname,
            source=cds_row.source,
            feature="exon",
            start=cds_row.start,
            end=cds_row.end,
            score=cds_row.score,
            strand=cds_row.strand,
            frame=".",
            attributes=exon_attris,
        )
        if strand == "-" and i == 0:
            exon_row.start -= 3
        elif strand != "-" and i == len(cds_rows) - 1:
            exon_row.end += 3
        exon_rows.append(exon_row)
    return exon_rows


def merge_overlapping_rows(rows, kind):
    tid = rows[0].attributes["transcript_id"]
    i = 0
    while i < len(rows) - 1:
        row1, row2 = rows[i], rows[i + 1]
        if row1.end + 1 >= row2.start:
            LOGGER.warning(
                f"Transcript {tid} has overlapping {kind}s: "
                f"{row1.start}-{row1.end} and {row2.start}-{row2.end}. Merging."
            )
            row1.end = max(row1.end, row2.end)
            rows.pop(i + 1)
        else:
            i += 1
    return rows


def _cds_exon_offset(exon_rows, cds_rows):
    """Return the exon index containing the first CDS, or None."""
    for ei, exon_row in enumerate(exon_rows):
        if exon_row.end < cds_rows[0].start:
            continue
        if exon_row.start > cds_rows[0].end:
            return None
        return ei
    return None


def check_cds_consistency_with_exon(exon_rows, cds_rows):
    offset = _cds_exon_offset(exon_rows, cds_rows)
    if offset is None or offset + len(cds_rows) > len(exon_rows):
        return False
    for ci, cds_row in enumerate(cds_rows):
        exon_row = exon_rows[offset + ci]
        if ci == 0:
            if ci == len(cds_rows) - 1:
                if not (
                    cds_row.start >= exon_row.start and cds_row.end <= exon_row.end
                ):
                    return False
            elif not (cds_row.start >= exon_row.start and cds_row.end == exon_row.end):
                return False
        else:
            if ci == len(cds_rows) - 1:
                if not (
                    cds_row.start == exon_row.start and cds_row.end <= exon_row.end
                ):
                    return False
            elif not (cds_row.start == exon_row.start and cds_row.end == exon_row.end):
                return False
    return True


def assign_exon_number_for_exon_rows(exon_rows):
    strand = exon_rows[0].strand
    for i, exon_row in enumerate(exon_rows):
        exon_row.attributes["exon_number"] = (
            (len(exon_rows) - i) if strand == "-" else (i + 1)
        )


def assign_exon_number_for_cds_rows(exon_rows, cds_rows):
    offset = _cds_exon_offset(exon_rows, cds_rows)
    if offset is None:
        return False
    for ci, cds_row in enumerate(cds_rows):
        cds_row.attributes["exon_number"] = exon_rows[offset + ci].attributes[
            "exon_number"
        ]
    return True


def get_canonical_transcript_id(final_gene_rows, canonicals=None):
    transcript_lengths = defaultdict(int)
    for row in final_gene_rows:
        if row.feature == "exon":
            transcript_lengths[row.attributes["transcript_id"]] += row.length

    tids1, tids2, tids3, tids4 = [], [], [], []
    for row in final_gene_rows:
        if row.feature != "transcript":
            continue
        transcript_id = row.attributes["transcript_id"]
        transcript_type = row.attributes["transcript_type"]
        if canonicals is not None and (
            transcript_id in canonicals or transcript_id.split(".")[0] in canonicals
        ):
            (tids1 if transcript_type in ("protein_coding", "mRNA") else tids2).append(
                transcript_id
            )
        (tids3 if transcript_type in ("protein_coding", "mRNA") else tids4).append(
            transcript_id
        )
    vs = tids1 or tids2 or tids3 or tids4
    if not vs:
        return None
    return sorted(vs, key=lambda x: transcript_lengths[x], reverse=True)[0]


def _add_ensembl_canonical_tag(attributes):
    """Mark a transcript row so `prepare`'s ranking puts it first."""
    tags = [t for t in attributes.get("tag", "").split("; ") if t]
    if "Ensembl_canonical" not in tags:
        tags.append("Ensembl_canonical")
        attributes["tag"] = "; ".join(tags)


# ---------------------------------------------------------------------------
# FASTA refiner
# ---------------------------------------------------------------------------


class FastaRefiner:
    def __init__(
        self, input_fasta, output_prefix, seqname_mapper=None, seqname_pattern=None
    ):
        self.input_fasta = input_fasta
        self.output_fasta = output_prefix + ".genome.fasta"
        self.output_sizes = output_prefix + ".genome.sizes"
        self.seqname_mapper = seqname_mapper
        self.seqname_pattern = seqname_pattern

    def run(self):
        LOGGER.info(f"Refining genome fasta: {self.input_fasta}")
        f = (
            gzip.open(self.input_fasta, "rt")
            if self.input_fasta.endswith(".gz")
            else open(self.input_fasta, "r")
        )
        n_kept = 0
        with open(self.output_fasta, "w") as fw:
            keep = False
            for line in f:
                if line.startswith(">"):
                    s = line.strip("\n")[1:].strip()
                    i = s.find(" ")
                    if i == -1:
                        seqname_old, description = s, None
                    else:
                        seqname_old, description = s[:i], s[i + 1 :].strip()

                    seqname_new = seqname_old
                    if self.seqname_mapper is not None:
                        if seqname_old in self.seqname_mapper:
                            seqname_new = self.seqname_mapper[seqname_old]
                        else:
                            LOGGER.warning(
                                f"{seqname_old} is not in the seqname mapper."
                            )

                    if self.seqname_pattern is not None:
                        if re.search(self.seqname_pattern, seqname_new) is None:
                            keep = False
                            continue

                    line = (
                        f">{seqname_new} {seqname_old}\n"
                        if description is None
                        else f">{seqname_new} {seqname_old} {description}\n"
                    )
                    fw.write(line)
                    keep = True
                    n_kept += 1
                elif keep:
                    fw.write(line)
        f.close()
        LOGGER.info(f"Kept {n_kept} sequences")

        # faidx + sizes via pysam (no external samtools needed)
        pysam.faidx(self.output_fasta)
        with (
            open(self.output_fasta + ".fai") as fai,
            open(self.output_sizes, "w") as fw,
        ):
            for line in fai:
                row = line.split("\t")
                fw.write(f"{row[0]}\t{row[1]}\n")


# ---------------------------------------------------------------------------
# GTF refiner
# ---------------------------------------------------------------------------


class GtfRefiner:
    def __init__(
        self,
        input_gtf,
        output_prefix,
        seqname_mapper=None,
        seqname_pattern=None,
        canonicals=None,
    ):
        self.input_gtf = input_gtf
        self.output_gtf = output_prefix + ".annotation.gtf"
        self.output_skip_gtf = output_prefix + ".skip.gtf"
        self.feature_summary_txt = output_prefix + ".gene_features_summary.txt"
        self.seqname_mapper = seqname_mapper
        self.seqname_pattern = seqname_pattern
        self.canonicals = canonicals

        self.exist_gene_names = set()
        self.exist_transcript_ids = set()
        self.exist_transcript_names = set()

    # -- driver ----------------------------------------------------------------

    def run(self):
        LOGGER.info(f"Loading annotation from {self.input_gtf}")
        annotation = GeneAnnotation.from_file(
            self.input_gtf,
            seqname_mapper=self.seqname_mapper,
            seqname_pattern=self.seqname_pattern,
        )
        dropped = annotation.prune()
        kept = sum(len(g.rows) for g in annotation.iter_genes())
        LOGGER.info(
            f"Loaded {kept} rows ({dropped} dropped, {len(annotation.genes)} genes)"
        )

        total_genes = len(annotation.genes)
        features_counter = defaultdict(int)
        succeeded = 0
        skipped_rows = []
        for gene in list(annotation.iter_genes()):
            features_counter[tuple(sorted({r.feature for r in gene.rows}))] += 1
            if self.process_gene(gene):
                succeeded += 1
            else:
                skipped_rows.extend(list(gene.rows))
                annotation.remove_gene(gene.gene_id)

        if total_genes:
            LOGGER.info(
                f"Refined {succeeded}/{total_genes} genes "
                f"({succeeded * 100.0 / total_genes:.2f}%)"
            )
        else:
            LOGGER.warning("No usable genes found in the input GTF.")

        annotation.write_gtf(
            self.output_gtf,
            extra_comments=("!refined gtf",),
            sort=True,
            check=True,
            bgzip=True,
        )
        with open(self.output_skip_gtf, "w") as fw:
            for row in skipped_rows:
                fw.write(row.to_gtf_line() + "\n")

        total = sum(features_counter.values())
        with open(self.feature_summary_txt, "w") as fw:
            fw.write("Total\tRatio\tFeatures\n")
            if total:
                for features, count in features_counter.items():
                    fw.write(f"{count}\t{count / total:.6f}\t{','.join(features)}\n")

    # -- per-gene transform -----------------------------------------------------

    def process_gene(self, gene):
        rows = gene.rows
        gene_id = gene.gene_id

        if not check_seqname_and_strand_consistency(rows):
            LOGGER.warning(f"Gene {gene_id} has inconsistent rows. Skipping.")
            return False

        gene_type = infer_gene_type(rows)
        gene_name = infer_gene_name(rows)

        if gene_name in self.exist_gene_names:
            gene_name = f"{gene_name}_{gene_id}"
            LOGGER.warning(
                f"Gene {gene_id} has a duplicate gene name; renamed to {gene_name}"
            )
        self.exist_gene_names.add(gene_name)

        for row in rows:
            row.attributes["gene_type"] = gene_type
            # biotype aliases: `prepare` reads the Ensembl names; GENCODE input
            # only carries gene_type, so keep both in sync for compatibility.
            row.attributes["gene_biotype"] = gene_type
            row.attributes["gene_name"] = gene_name

        data_features = group_rows_by_feature(rows)
        if "gene" in data_features:
            if len(data_features["gene"]) > 1:
                LOGGER.warning(f"Gene {gene_id} has multiple gene features. Skipping.")
                return False
            gene_row = data_features["gene"][0]
        else:
            LOGGER.warning(f"Gene {gene_id} has no gene feature; creating one.")
            gene_row = create_gene_row_from_gene_rows(rows)
        if gene_row is None:
            LOGGER.warning(f"Gene {gene_id}: cannot create a gene row. Skipping.")
            return False
        final_gene_rows = [gene_row]

        if not gene.transcripts:
            LOGGER.warning(
                f"Gene {gene_id} has no transcript; creating transcript + exon."
            )
            transcript_row = create_transcript_row_from_gene_row(gene_row)
            final_gene_rows.append(transcript_row)
            final_gene_rows.append(create_exon_row_from_transcript_row(transcript_row))
        else:
            for tid in list(gene.transcripts):
                final_transcript_rows = self.process_transcript(gene.transcripts[tid])
                if final_transcript_rows is None:
                    return False
                final_gene_rows.extend(final_transcript_rows)

        canonical_transcript_id = get_canonical_transcript_id(
            final_gene_rows, self.canonicals
        )
        for row in final_gene_rows:
            if row.feature == "gene":
                continue
            row.attributes["is_canonical"] = (
                row.attributes["transcript_id"] == canonical_transcript_id
            )
            if (
                canonical_transcript_id
                and row.attributes["transcript_id"] == canonical_transcript_id
            ):
                _add_ensembl_canonical_tag(row.attributes)

        gene.rows = final_gene_rows
        return True

    def process_transcript(self, tx):
        rows = tx.rows
        tid = tx.transcript_id

        if tid in self.exist_transcript_ids:
            tid = f"{tid}_{tx.gene_id}"
            LOGGER.warning(f"Duplicate transcript id; renamed to {tid}")
            for row in rows:
                row.attributes["transcript_id"] = tid
        self.exist_transcript_ids.add(tid)

        transcript_type = infer_transcript_type(rows, tid)
        transcript_name = infer_transcript_name(rows, tid)
        if transcript_name in self.exist_transcript_names:
            transcript_name = f"{transcript_name}_{tid}"
            LOGGER.warning(f"Duplicate transcript name; renamed to {transcript_name}")
        self.exist_transcript_names.add(transcript_name)

        for row in rows:
            row.attributes["transcript_type"] = transcript_type
            # biotype alias for `prepare` (see process_gene)
            row.attributes["transcript_biotype"] = transcript_type
            row.attributes["transcript_name"] = transcript_name

        data_features = group_rows_by_feature(rows)

        if "transcript" in data_features:
            if len(data_features["transcript"]) > 1:
                LOGGER.error(f"Transcript {tid} has multiple transcript features.")
                return None
            transcript_row = data_features["transcript"][0]
        else:
            transcript_row = create_transcript_row_from_transcript_rows(rows)
        if transcript_row is None:
            LOGGER.error(f"Transcript {tid}: cannot create a transcript row.")
            return None
        final_transcript_rows = [transcript_row]

        for row in rows:
            if not (
                transcript_row.start <= row.start and row.end <= transcript_row.end
            ):
                LOGGER.error(
                    f"Transcript {tid} has rows outside the transcript bounds."
                )
                return None

        cds_rows = None
        if "CDS" in data_features:
            cds_rows = sorted(data_features["CDS"], key=lambda r: (r.start, r.end))
            cds_rows = merge_overlapping_rows(cds_rows, "CDS")

        if "exon" in data_features:
            exon_rows = sorted(data_features["exon"], key=lambda r: (r.start, r.end))
            exon_rows = merge_overlapping_rows(exon_rows, "exon")
        elif "CDS" in data_features:
            exon_rows = create_exon_rows_from_cds_rows(cds_rows)
        else:
            exon_rows = [create_exon_row_from_transcript_row(transcript_row)]
        assign_exon_number_for_exon_rows(exon_rows)
        final_transcript_rows.extend(exon_rows)

        if "CDS" in data_features:
            if not check_cds_consistency_with_exon(exon_rows, cds_rows):
                LOGGER.error(f"Transcript {tid}: CDS rows do not match exons.")
                return None
            assign_exon_number_for_cds_rows(exon_rows, cds_rows)
            final_transcript_rows.extend(cds_rows)

        # carry codon/UTR/other feature rows through (drop their exon_number;
        # codon/UTR rows are re-numbered by the exons when present)
        for feature in data_features:
            if feature in ("gene", "transcript", "exon", "CDS"):
                continue
            for row in data_features[feature]:
                row.attributes.pop("exon_number", None)
                final_transcript_rows.append(row)

        return final_transcript_rows


# ---------------------------------------------------------------------------
# public entry point
# ---------------------------------------------------------------------------


def refine_genome_references(
    input_fasta=None,
    input_gtf=None,
    outdir="./",
    name=None,
    rename_mapper=None,
    seqname_pattern=None,
    canonical_transcripts=None,
):
    """Refine a genome FASTA and/or GTF for downstream coralsnake commands."""
    if input_fasta is None and input_gtf is None:
        raise ValueError("Nothing to refine: pass --fasta-file and/or --gtf-file.")

    if name is None:
        name = os.path.basename(os.path.abspath(outdir))
    if not name:
        name = "refined"
    os.makedirs(outdir, exist_ok=True)
    output_prefix = os.path.join(outdir, name)

    seqname_mapper = (
        None if rename_mapper is None else load_seqname_mapper(rename_mapper)
    )
    canonicals = (
        None
        if canonical_transcripts is None
        else load_canonical_transcripts(canonical_transcripts)
    )

    if input_fasta is not None:
        FastaRefiner(
            input_fasta=input_fasta,
            output_prefix=output_prefix,
            seqname_mapper=seqname_mapper,
            seqname_pattern=seqname_pattern,
        ).run()

    if input_gtf is not None:
        GtfRefiner(
            input_gtf=input_gtf,
            output_prefix=output_prefix,
            seqname_mapper=seqname_mapper,
            seqname_pattern=seqname_pattern,
            canonicals=canonicals,
        ).run()

    return output_prefix
