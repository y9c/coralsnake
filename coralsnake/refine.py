#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Reference cleaning for coralsnake (`coralsnake refine`).
#
# Rewritten on the package infrastructure:
#   * attribute parsing reuses gtf2tx.parse_gtf_annot (regex-based, handles
#     quoted ';' inside values and missing trailing ';')
#   * faidx / bgzip / tabix use pysam (no external samtools/bgzip/tabix on PATH)
#   * the refined GTF keeps start/stop codon, UTR and other feature rows
#   * confirmed crashes fixed (exon-only genes, CDS without a containing exon,
#     empty default name)

import gzip
import os
import re
from collections import defaultdict

import pysam

from .gtf2tx import parse_gtf_annot
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


class FastaRefiner:
    def __init__(self, input_fasta, output_prefix, seqname_mapper=None, seqname_pattern=None):
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
                            LOGGER.warning(f"{seqname_old} is not in the seqname mapper.")

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
        with open(self.output_fasta + ".fai") as fai, open(self.output_sizes, "w") as fw:
            for line in fai:
                row = line.split("\t")
                fw.write(f"{row[0]}\t{row[1]}\n")


class GtfUtils:
    @staticmethod
    def group_rows_by_gene_id(rows):
        data = defaultdict(list)
        for row in rows:
            data[row[8]["gene_id"]].append(row)
        return data

    @staticmethod
    def group_rows_by_transcript_id(rows):
        data = defaultdict(list)
        for row in rows:
            if row[2] == "gene":
                continue
            data[row[8]["transcript_id"]].append(row)
        return data

    @staticmethod
    def group_rows_by_feature(rows):
        data = defaultdict(list)
        for row in rows:
            data[row[2]].append(row)
        return data

    @staticmethod
    def check_seqname_and_strand_consistency(rows):
        row1 = rows[0]
        gid = row1[8]["gene_id"]
        for row2 in rows[1:]:
            if row1[0] != row2[0]:
                LOGGER.error(
                    f"Gene {gid} has rows on different seqnames: {row1[0]} and {row2[0]}"
                )
                return False
            if row1[6] != row2[6]:
                LOGGER.error(
                    f"Gene {gid} has rows on different strands: {row1[6]} and {row2[6]}"
                )
                return False
        return True

    @staticmethod
    def load_gtf(path, seqname_mapper=None, seqname_pattern=None):
        header = []
        rows = []
        f = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")
        for line in f:
            if line.startswith("#"):
                header.append(line.strip())
                continue
            row = line.strip().split("\t")
            if len(row) != 9:
                LOGGER.warning(f"Skipping malformed GTF line ({len(row)} fields): {line.strip()[:80]}")
                continue
            if seqname_mapper is not None:
                if row[0] in seqname_mapper:
                    row[0] = seqname_mapper[row[0]]
                else:
                    LOGGER.warning(f"{row[0]} is not in the seqname mapper.")
            if seqname_pattern is not None and re.search(seqname_pattern, row[0]) is None:
                continue
            row[1] = row[1].replace(" ", "_")
            row[3] = int(row[3])
            row[4] = int(row[4])
            if row[3] < 1 or row[4] < row[3]:
                LOGGER.warning(f"Skipping line with invalid coordinates: {line.strip()[:80]}")
                continue
            row[8] = parse_gtf_annot(row[8])
            if "gene_id" not in row[8]:
                LOGGER.warning(f"Skipping line without gene_id: {line.strip()[:80]}")
                continue
            if row[2] == "gene" and "transcript_id" in row[8]:
                del row[8]["transcript_id"]
            if row[2] != "gene" and "transcript_id" not in row[8]:
                LOGGER.warning(f"Skipping non-gene line without transcript_id: {line.strip()[:80]}")
                continue
            rows.append(row)
        f.close()
        return header, rows

    @staticmethod
    def _bounding_row(rows, feature, drop_transcript_id=True):
        """Build a synthetic row spanning all ``rows`` (min start, max end)."""
        gene_row = None
        attris = None
        for row in rows:
            if attris is None:
                attris = row[8].copy()
                if drop_transcript_id:
                    attris.pop("transcript_id", None)
            if gene_row is None:
                gene_row = row[:8]
                gene_row[2] = feature
                gene_row.append(attris)
            else:
                gene_row[3] = min(gene_row[3], row[3])
                gene_row[4] = max(gene_row[4], row[4])
        return gene_row

    @staticmethod
    def create_gene_row_from_gene_rows(rows):
        gene_id = rows[0][8]["gene_id"]
        data_features = GtfUtils.group_rows_by_feature(rows)
        if "transcript" in data_features:
            LOGGER.info(f"Creating gene row from transcript rows for {gene_id}")
            return GtfUtils._bounding_row(data_features["transcript"], "gene")
        if "exon" in data_features:
            LOGGER.info(f"Creating gene row from exon rows for {gene_id}")
            return GtfUtils._bounding_row(data_features["exon"], "gene")
        if "CDS" in data_features:
            LOGGER.info(f"Creating gene row from CDS rows for {gene_id}")
            gene_row = GtfUtils._bounding_row(data_features["CDS"], "gene")
            gene_row[7] = "."
            # extend to include the stop codon
            if gene_row[6] == "-":
                gene_row[3] -= 3
            else:
                gene_row[4] += 3
            return gene_row
        return None

    @staticmethod
    def _infer_unique(rows, gid, keys, kind):
        """Collect a unique attribute value across rows; None when ambiguous."""
        vs = set()
        for row in rows:
            attris = row[8]
            for k in keys:
                if k in attris:
                    vs.add(attris[k])
        if len(vs) == 1:
            return vs.pop()
        if len(vs) > 1:
            LOGGER.info(f"{kind} {gid} has multiple values for {keys}: {sorted(vs)}")
        return None

    @staticmethod
    def infer_gene_type(rows):
        gid = rows[0][8]["gene_id"]
        v = GtfUtils._infer_unique(rows, gid, ["gene_type", "gene_biotype"], "Gene")
        if v is not None:
            return v
        for row in rows:
            if row[2] == "CDS":
                return "protein_coding"
        return "unknown"

    @staticmethod
    def infer_gene_name(rows):
        gid = rows[0][8]["gene_id"]
        v = GtfUtils._infer_unique(rows, gid, ["gene_name", "gene"], "Gene")
        return v if v is not None else gid

    @staticmethod
    def infer_transcript_type(rows):
        tid = rows[0][8]["transcript_id"]
        v = GtfUtils._infer_unique(rows, tid, ["transcript_type", "transcript_biotype"], "Transcript")
        if v is not None:
            return v
        for row in rows:
            if row[2] == "CDS":
                return "protein_coding"
        return "unknown"

    @staticmethod
    def infer_transcript_name(rows):
        tid = rows[0][8]["transcript_id"]
        v = GtfUtils._infer_unique(rows, tid, ["transcript_name"], "Transcript")
        return v if v is not None else tid

    @staticmethod
    def output_row(fw, row, check=False):
        attris = row[8]
        if check:
            assert "gene_id" in attris
            assert "gene_name" in attris
            assert "gene_type" in attris
            if row[2] == "gene":
                assert "transcript_id" not in attris
            else:
                assert "transcript_id" in attris
                assert "transcript_name" in attris
                assert "transcript_type" in attris
        items = [f'{k} "{v}";' for k, v in attris.items()]
        fw.write("\t".join(str(x) for x in row[:8]) + "\t" + " ".join(items) + "\n")

    @staticmethod
    def output_rows(fw, rows, check=False):
        for row in rows:
            GtfUtils.output_row(fw, row, check)

    @staticmethod
    def create_transcript_row_from_transcript_rows(rows):
        data_features = GtfUtils.group_rows_by_feature(rows)
        if "exon" in data_features:
            return GtfUtils._bounding_row(data_features["exon"], "transcript", drop_transcript_id=False)
        if "CDS" in data_features:
            transcript_row = GtfUtils._bounding_row(data_features["CDS"], "transcript", drop_transcript_id=False)
            transcript_row[7] = "."
            if transcript_row[6] == "-":
                transcript_row[3] -= 3
            else:
                transcript_row[4] += 3
            return transcript_row
        return None

    @staticmethod
    def create_transcript_row_from_gene_row(gene_row):
        gene_attris = gene_row[8]
        transcript_id = gene_attris["gene_id"] + "_transcript"
        transcript_row = gene_row[:8]
        transcript_row[2] = "transcript"
        transcript_attris = gene_attris.copy()
        transcript_attris["transcript_id"] = transcript_id
        transcript_attris["transcript_name"] = transcript_id
        transcript_attris["transcript_type"] = gene_attris["gene_type"]
        transcript_row.append(transcript_attris)
        return transcript_row

    @staticmethod
    def create_exon_row_from_transcript_row(transcript_row):
        exon_row = transcript_row[:8]
        exon_row[2] = "exon"
        exon_attris = transcript_row[8].copy()
        exon_attris["exon_number"] = 1
        exon_row.append(exon_attris)
        return exon_row

    @staticmethod
    def create_exon_rows_from_cds_rows(cds_rows):
        exon_rows = []
        strand = cds_rows[0][6]
        for i, cds_row in enumerate(cds_rows):
            exon_row = cds_row[:8]
            exon_row[2] = "exon"
            exon_row[7] = "."
            exon_row.append(cds_row[8].copy())
            if strand == "-" and i == 0:
                exon_row[3] -= 3
            elif strand != "-" and i == len(cds_rows) - 1:
                exon_row[4] += 3
            exon_rows.append(exon_row)
        return exon_rows

    @staticmethod
    def merge_overlapping_rows(rows, kind):
        tid = rows[0][8]["transcript_id"]
        i = 0
        while i < len(rows) - 1:
            row1, row2 = rows[i], rows[i + 1]
            if row1[4] + 1 >= row2[3]:
                LOGGER.warning(
                    f"Transcript {tid} has overlapping {kind}s: "
                    f"{row1[3]}-{row1[4]} and {row2[3]}-{row2[4]}. Merging."
                )
                row1[4] = max(row1[4], row2[4])
                rows.pop(i + 1)
            else:
                i += 1
        return rows

    @staticmethod
    def _cds_exon_offset(exon_rows, cds_rows):
        """Return the exon index containing the first CDS, or None."""
        for ei, exon_row in enumerate(exon_rows):
            if exon_row[4] < cds_rows[0][3]:
                continue
            if exon_row[3] > cds_rows[0][4]:
                return None
            return ei
        return None

    @staticmethod
    def check_cds_consistency_with_exon(exon_rows, cds_rows):
        offset = GtfUtils._cds_exon_offset(exon_rows, cds_rows)
        if offset is None or offset + len(cds_rows) > len(exon_rows):
            return False
        for ci, cds_row in enumerate(cds_rows):
            exon_row = exon_rows[offset + ci]
            if ci == 0:
                if ci == len(cds_rows) - 1:
                    if not (cds_row[3] >= exon_row[3] and cds_row[4] <= exon_row[4]):
                        return False
                elif not (cds_row[3] >= exon_row[3] and cds_row[4] == exon_row[4]):
                    return False
            else:
                if ci == len(cds_rows) - 1:
                    if not (cds_row[3] == exon_row[3] and cds_row[4] <= exon_row[4]):
                        return False
                elif not (cds_row[3] == exon_row[3] and cds_row[4] == exon_row[4]):
                    return False
        return True

    @staticmethod
    def assign_exon_number_for_exon_rows(exon_rows):
        strand = exon_rows[0][6]
        for i, exon_row in enumerate(exon_rows):
            exon_row[8]["exon_number"] = (len(exon_rows) - i) if strand == "-" else (i + 1)

    @staticmethod
    def assign_exon_number_for_cds_rows(exon_rows, cds_rows):
        offset = GtfUtils._cds_exon_offset(exon_rows, cds_rows)
        if offset is None:
            return False
        for ci, cds_row in enumerate(cds_rows):
            cds_row[8]["exon_number"] = exon_rows[offset + ci][8]["exon_number"]
        return True

    # Feature rows carried through unchanged (besides attribute normalization)
    @staticmethod
    def get_canonical_transcript_id(final_gene_rows, canonicals=None):
        transcript_lengths = defaultdict(int)
        for row in final_gene_rows:
            if row[2] == "exon":
                transcript_lengths[row[8]["transcript_id"]] += row[4] - row[3] + 1

        tids1, tids2, tids3, tids4 = [], [], [], []
        for row in final_gene_rows:
            if row[2] != "transcript":
                continue
            transcript_id = row[8]["transcript_id"]
            transcript_type = row[8]["transcript_type"]
            if canonicals is not None and (
                transcript_id in canonicals or transcript_id.split(".")[0] in canonicals
            ):
                (tids1 if transcript_type in ("protein_coding", "mRNA") else tids2).append(transcript_id)
            (tids3 if transcript_type in ("protein_coding", "mRNA") else tids4).append(transcript_id)
        vs = tids1 or tids2 or tids3 or tids4
        if not vs:
            return None
        return sorted(vs, key=lambda x: transcript_lengths[x], reverse=True)[0]


class GtfRefiner:
    def __init__(self, input_gtf, output_prefix, seqname_mapper=None, seqname_pattern=None, canonicals=None):
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

    def run(self):
        with open(self.output_gtf, "w") as self.fw, open(self.output_skip_gtf, "w") as self.fw_skip:
            LOGGER.info(f"Loading rows from {self.input_gtf}")
            header, rows = GtfUtils.load_gtf(
                self.input_gtf,
                seqname_mapper=self.seqname_mapper,
                seqname_pattern=self.seqname_pattern,
            )
            LOGGER.info(f"Loaded {len(rows)} rows")
            for line in header:
                self.fw.write(line + "\n")
            self.fw.write("#!refined gtf\n")

            data_gene = GtfUtils.group_rows_by_gene_id(rows)
            total_genes = len(data_gene)
            LOGGER.info(f"Grouped rows into {total_genes} genes")

            features_counter = defaultdict(int)
            succeed_genes = 0
            for gene_rows in data_gene.values():
                features_counter[tuple(sorted({r[2] for r in gene_rows}))] += 1
                if self.process_gene_rows(gene_rows):
                    succeed_genes += 1
            LOGGER.info(
                f"Refined {succeed_genes}/{total_genes} genes "
                f"({succeed_genes * 100.0 / total_genes:.2f}%)"
            )

        self._sort_and_index()
        with open(self.feature_summary_txt, "w") as fw:
            total = sum(features_counter.values())
            fw.write("Total\tRatio\tFeatures\n")
            for features, count in features_counter.items():
                fw.write(f"{count}\t{count / total:.6f}\t{','.join(features)}\n")

    def _sort_and_index(self):
        """Sort the GTF by (seqname, start, end) and bgzip+tabix via pysam."""
        with open(self.output_gtf) as f:
            header = [ln for ln in f if ln.startswith("#")]
        with open(self.output_gtf) as f:
            body = [ln for ln in f if not ln.startswith("#")]

        def _key(line):
            row = line.split("\t")
            return (row[0], int(row[3]), int(row[4]))

        body.sort(key=_key)
        sorted_gtf = self.output_gtf
        with open(sorted_gtf, "w") as f:
            f.writelines(header + body)
        try:
            gz = sorted_gtf + ".gz"
            pysam.tabix_compress(sorted_gtf, gz, force=True)
            pysam.tabix_index(gz, preset="gff", force=True)
        except Exception as e:  # pragma: no cover - tabix is best-effort
            LOGGER.warning(f"Skipping bgzip/tabix indexing: {e}")

    def process_gene_rows(self, rows):
        final_gene_rows = []
        gene_id = rows[0][8]["gene_id"]

        if not GtfUtils.check_seqname_and_strand_consistency(rows):
            LOGGER.warning(f"Gene {gene_id} has inconsistent rows. Skipping.")
            GtfUtils.output_rows(self.fw_skip, rows)
            return False

        gene_type = GtfUtils.infer_gene_type(rows)
        gene_name = GtfUtils.infer_gene_name(rows)

        if gene_name in self.exist_gene_names:
            gene_name = f"{gene_name}_{gene_id}"
            LOGGER.warning(f"Gene {gene_id} has a duplicate gene name; renamed to {gene_name}")
        self.exist_gene_names.add(gene_name)

        for row in rows:
            row[8]["gene_type"] = gene_type
            row[8]["gene_name"] = gene_name

        data_features = GtfUtils.group_rows_by_feature(rows)
        if "gene" in data_features:
            if len(data_features["gene"]) > 1:
                LOGGER.warning(f"Gene {gene_id} has multiple gene features. Skipping.")
                GtfUtils.output_rows(self.fw_skip, rows)
                return False
            gene_row = data_features["gene"][0]
        else:
            LOGGER.warning(f"Gene {gene_id} has no gene feature; creating one.")
            gene_row = GtfUtils.create_gene_row_from_gene_rows(rows)
        if gene_row is None:
            LOGGER.warning(f"Gene {gene_id}: cannot create a gene row. Skipping.")
            GtfUtils.output_rows(self.fw_skip, rows)
            return False
        final_gene_rows.append(gene_row)

        data_transcript = GtfUtils.group_rows_by_transcript_id(rows)
        if len(data_transcript) == 0:
            LOGGER.warning(f"Gene {gene_id} has no transcript; creating transcript + exon.")
            transcript_row = GtfUtils.create_transcript_row_from_gene_row(gene_row)
            final_gene_rows.append(transcript_row)
            final_gene_rows.append(GtfUtils.create_exon_row_from_transcript_row(transcript_row))
        else:
            for tx_rows in data_transcript.values():
                final_transcript_rows = self.process_transcript_rows(tx_rows)
                if final_transcript_rows is None:
                    GtfUtils.output_rows(self.fw_skip, rows)
                    return False
                final_gene_rows.extend(final_transcript_rows)

        canonical_transcript_id = GtfUtils.get_canonical_transcript_id(final_gene_rows, self.canonicals)
        for row in final_gene_rows:
            if row[2] != "gene":
                row[8]["is_canonical"] = row[8]["transcript_id"] == canonical_transcript_id
            GtfUtils.output_row(self.fw, row, check=True)
        return True

    def process_transcript_rows(self, rows):
        final_transcript_rows = []
        tid = rows[0][8]["transcript_id"]

        if tid in self.exist_transcript_ids:
            gid = rows[0][8]["gene_id"]
            tid = f"{tid}_{gid}"
            LOGGER.warning(f"Duplicate transcript id; renamed to {tid}")
            for row in rows:
                row[8]["transcript_id"] = tid
        self.exist_transcript_ids.add(tid)

        transcript_type = GtfUtils.infer_transcript_type(rows)
        transcript_name = GtfUtils.infer_transcript_name(rows)
        if transcript_name in self.exist_transcript_names:
            transcript_name = f"{transcript_name}_{tid}"
            LOGGER.warning(f"Duplicate transcript name; renamed to {transcript_name}")
        self.exist_transcript_names.add(transcript_name)

        for row in rows:
            row[8]["transcript_type"] = transcript_type
            row[8]["transcript_name"] = transcript_name

        data_features = GtfUtils.group_rows_by_feature(rows)

        if "transcript" in data_features:
            if len(data_features["transcript"]) > 1:
                LOGGER.error(f"Transcript {tid} has multiple transcript features.")
                return None
            transcript_row = data_features["transcript"][0]
        else:
            transcript_row = GtfUtils.create_transcript_row_from_transcript_rows(rows)
        if transcript_row is None:
            LOGGER.error(f"Transcript {tid}: cannot create a transcript row.")
            return None
        final_transcript_rows.append(transcript_row)

        for row in rows:
            if not (transcript_row[3] <= row[3] and row[4] <= transcript_row[4]):
                LOGGER.error(f"Transcript {tid} has rows outside the transcript bounds.")
                return None

        cds_rows = None
        if "CDS" in data_features:
            cds_rows = sorted(data_features["CDS"], key=lambda row: [row[3], row[4]])
            cds_rows = GtfUtils.merge_overlapping_rows(cds_rows, "CDS")

        if "exon" in data_features:
            exon_rows = sorted(data_features["exon"], key=lambda row: [row[3], row[4]])
            exon_rows = GtfUtils.merge_overlapping_rows(exon_rows, "exon")
        elif "CDS" in data_features:
            exon_rows = GtfUtils.create_exon_rows_from_cds_rows(cds_rows)
        else:
            exon_rows = [GtfUtils.create_exon_row_from_transcript_row(transcript_row)]
        GtfUtils.assign_exon_number_for_exon_rows(exon_rows)
        final_transcript_rows.extend(exon_rows)

        if "CDS" in data_features:
            if not GtfUtils.check_cds_consistency_with_exon(exon_rows, cds_rows):
                LOGGER.error(f"Transcript {tid}: CDS rows do not match exons.")
                return None
            GtfUtils.assign_exon_number_for_cds_rows(exon_rows, cds_rows)
            final_transcript_rows.extend(cds_rows)

        # carry codon/UTR/other feature rows through (drop their exon_number;
        # codon/UTR rows are re-numbered by the exons when present)
        for feature in data_features:
            if feature in ("gene", "transcript", "exon", "CDS"):
                continue
            for row in data_features[feature]:
                row[8].pop("exon_number", None)
                final_transcript_rows.append(row)

        return final_transcript_rows


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

    seqname_mapper = None if rename_mapper is None else load_seqname_mapper(rename_mapper)
    canonicals = (
        None if canonical_transcripts is None else load_canonical_transcripts(canonical_transcripts)
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
