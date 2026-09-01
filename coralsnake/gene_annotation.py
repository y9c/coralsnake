#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shared gene-annotation model for coralsnake.

One :class:`GeneAnnotation` object models one annotated genome: every GTF (or
GFF3-style) feature row is parsed once, kept losslessly, and grouped under its
gene and transcript. Coordinates stay in GTF convention (1-based, closed) and
0-based half-open spans are exposed as helpers.

Why this object exists
----------------------
The package has several GTF consumers (``prepare``/``gtf2tx``, ``refine``,
``metagene``/``annotate``). Historically each parse-path diverged (different
attribute regexes, different rules for malformed lines, different name
aliasing). ``GeneAnnotation`` is the single source of truth for *reading*
(regex, comments, seqname renaming/filtering, grouping) and *writing*
(attribute serialization, sorting, bgzip/tabix) a gene annotation:

* ``refine`` loads a ``GeneAnnotation``, mutates rows/attributes, and
  re-serializes with :meth:`GeneAnnotation.write_gtf`.
* ``gtf2tx.read_gtf`` consumes :meth:`GeneAnnotation.iter_rows` so `prepare`
  shares the exact same read layer (its per-row ranking / biotype decisions are
  unchanged).

High-throughput consumers (``metagene``/``annotate``) keep their vectorized
Polars loader for speed — the *semantics* they rely on (attribute tokens, which
features are exon/codon rows, name aliasing) are defined here.

Memory note
-----------
A ~3M-line human GTF loads in a few seconds and holds every row as a lightweight
dataclass. That is fine for one-shot cleaning/preparation (the metagene table is
already cached as parquet); if a future consumer needs streaming, row parsing is
isolated in :meth:`GeneAnnotation._iter_file_lines`.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from typing import Iterator, Optional

import pysam

from .utils import get_logger

LOGGER = get_logger(__name__)


# ---------------------------------------------------------------------------
# Attribute parsing (canonical definitions, also re-exported by gtf2tx)
# ---------------------------------------------------------------------------

# GTF attributes: `key "value"; key "value";`
gtf_pattern = re.compile(r'(\w+)\s+"(.*?)"(?:;|$)')
# GFF3 attributes: `key=value;key=value;`
gff_pattern = re.compile(r"(\w+)=([^;]*)")


def _parse_annot(annot, pattern):
    """Parse an attribute field with ``pattern`` into a dict.

    Duplicate keys are joined with ``"; "`` (GTF tags repeat, e.g.
    ``tag "basic"; tag "MANE_Select"``).
    """
    matches = pattern.findall(annot.rstrip("\n").rstrip(";"))
    d = {}
    for k, v in matches:
        if k in d:
            d[k] = f"{d[k]}; {v}"
        else:
            d[k] = v
    return d


def parse_gtf_annot(annot):
    """Parse a GTF attribute field into a dict.

    Robust to quoted ``;`` inside values and a missing trailing ``;``.
    """
    return _parse_annot(annot, gtf_pattern)


def parse_gff_annot(annot):
    """Parse a GFF3 attribute field into a dict."""
    return _parse_annot(annot, gff_pattern)


def detect_gff(path: str) -> bool:
    """Heuristic: does the file look like GFF3 (attribute parser style)?"""
    p = path.lower()
    return p.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz"))


# ---------------------------------------------------------------------------
# Row / record model
# ---------------------------------------------------------------------------


@dataclass
class AnnotationRow:
    """One feature line of a gene annotation (GTF columns 1-9).

    Coordinates follow the GTF convention: 1-based, closed (inclusive) — the
    whole package stores GTF coordinates this way and converts to 0-based
    half-open only at the boundary (see :attr:`span_0`).
    """

    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str = "."
    strand: str = "."
    frame: str = "."
    attributes: dict = field(default_factory=dict)

    @property
    def span_0(self) -> tuple[int, int]:
        """0-based, half-open equivalent of this row's interval."""
        return (self.start - 1, self.end)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def to_gtf_line(self) -> str:
        """Serialize back to a GTF line (attributes in insertion order)."""
        attrs = " ".join(f'{k} "{v}";' for k, v in self.attributes.items())
        return "\t".join(
            [
                self.seqname,
                self.source,
                self.feature,
                str(self.start),
                str(self.end),
                self.score,
                self.strand,
                self.frame,
                attrs,
            ]
        )

    def __str__(self) -> str:  # friendly, compact form for log messages
        return (
            f"{self.seqname}:{self.feature}:{self.start}-{self.end}"
            f"({self.attributes.get('gene_id', '-')}/"
            f"{self.attributes.get('transcript_id', '-')})"
        )


@dataclass
class Transcript:
    """All feature rows belonging to one transcript."""

    gene_id: str
    transcript_id: str
    rows: list[AnnotationRow] = field(default_factory=list)

    @property
    def seqname(self) -> str:
        return self.rows[0].seqname if self.rows else ""

    @property
    def strand(self) -> str:
        return self.rows[0].strand if self.rows else "."

    @property
    def features(self) -> dict[str, list[AnnotationRow]]:
        out: dict[str, list[AnnotationRow]] = {}
        for row in self.rows:
            out.setdefault(row.feature, []).append(row)
        return out

    @property
    def exons(self) -> list[AnnotationRow]:
        return self.features.get("exon", [])

    @property
    def cds(self) -> list[AnnotationRow]:
        return self.features.get("CDS", [])

    @property
    def start_codons(self) -> list[AnnotationRow]:
        return self.features.get("start_codon", [])

    @property
    def stop_codons(self) -> list[AnnotationRow]:
        return self.features.get("stop_codon", [])

    @property
    def utrs(self) -> list[AnnotationRow]:
        out = []
        for row in self.rows:
            if "utr" in row.feature.lower():
                out.append(row)
        return out

    @property
    def span(self) -> tuple[int, int]:
        """Bounding box in 1-based closed coordinates."""
        if not self.rows:
            return (0, 0)
        return (
            min(r.start for r in self.rows),
            max(r.end for r in self.rows),
        )

    @property
    def is_coding(self) -> bool:
        return any(r.feature == "CDS" for r in self.rows)


@dataclass
class Gene:
    """All feature rows belonging to one gene, plus its transcripts."""

    gene_id: str
    rows: list[AnnotationRow] = field(default_factory=list)
    transcripts: dict[str, Transcript] = field(default_factory=dict)

    @property
    def seqname(self) -> str:
        return self.rows[0].seqname if self.rows else ""

    @property
    def strand(self) -> str:
        return self.rows[0].strand if self.rows else "."

    @property
    def features(self) -> dict[str, list[AnnotationRow]]:
        out: dict[str, list[AnnotationRow]] = {}
        for row in self.rows:
            out.setdefault(row.feature, []).append(row)
        return out

    @property
    def gene_rows(self) -> list[AnnotationRow]:
        return self.features.get("gene", [])

    @property
    def span(self) -> tuple[int, int]:
        """Bounding box in 1-based closed coordinates."""
        if not self.rows:
            return (0, 0)
        return (
            min(r.start for r in self.rows),
            max(r.end for r in self.rows),
        )


# ---------------------------------------------------------------------------
# The annotation object
# ---------------------------------------------------------------------------


class GeneAnnotation:
    """A parsed gene annotation (GTF / GFF3-style attributes).

    Rows are kept losslessly (every attribute survives) and grouped under
    :attr:`genes` (by ``gene_id``) and per-:class:`Transcript` (by
    ``transcript_id``). Rows without a ``gene_id`` are kept in
    :attr:`unassigned_rows` (they are still visible via :meth:`iter_rows`); call
    :meth:`prune` to drop them like the GTF-cleaning workflow expects.
    """

    def __init__(
        self,
        path: Optional[str] = None,
        *,
        is_gff: Optional[bool] = None,
        seqname_mapper: Optional[dict[str, str]] = None,
        seqname_pattern: Optional[str] = None,
    ):
        self.path = path
        self.is_gff = detect_gff(path) if is_gff is None else bool(is_gff)
        self.seqname_mapper = seqname_mapper
        self.seqname_pattern = seqname_pattern

        self.header: list[str] = []  # original "#..." comment lines
        self.genes: dict[str, Gene] = {}
        self.unassigned_rows: list[AnnotationRow] = []
        self._all_rows: list[AnnotationRow] = []

        if path is not None:
            self._load()

    @classmethod
    def from_file(
        cls,
        path: str,
        *,
        is_gff: Optional[bool] = None,
        seqname_mapper: Optional[dict[str, str]] = None,
        seqname_pattern: Optional[str] = None,
    ) -> "GeneAnnotation":
        return cls(
            path,
            is_gff=is_gff,
            seqname_mapper=seqname_mapper,
            seqname_pattern=seqname_pattern,
        )

    # -- loading -----------------------------------------------------------

    def _iter_file_lines(self):
        """Yield (lineno, raw_line) with transparent gzip support."""
        f = (
            gzip.open(self.path, "rt")
            if self.path.endswith(".gz")
            else open(self.path, "r")
        )
        with f:
            for lineno, line in enumerate(f, 1):
                yield lineno, line

    def _munge_seqname(self, seqname: str) -> str:
        if self.seqname_mapper is not None:
            if seqname in self.seqname_mapper:
                seqname = self.seqname_mapper[seqname]
            else:
                LOGGER.warning(f"{seqname} is not in the seqname mapper.")
        return seqname

    def _load(self) -> None:
        parse_annot = parse_gff_annot if self.is_gff else parse_gtf_annot
        rows: list[AnnotationRow] = []
        for lineno, line in self._iter_file_lines():
            if line.startswith("#"):
                stripped = line.strip()
                if stripped == "##FASTA" or stripped.startswith("##FASTA "):
                    break  # embedded sequence section; nothing else to parse
                self.header.append(stripped)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                LOGGER.warning(
                    f"Line {lineno}: skipping malformed annotation line "
                    f"({len(fields)} fields): {line.strip()[:80]}"
                )
                continue
            seqname = self._munge_seqname(fields[0])
            if (
                self.seqname_pattern is not None
                and re.search(self.seqname_pattern, seqname) is None
            ):
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                LOGGER.warning(
                    f"Line {lineno}: skipping line with non-integer coordinates: "
                    f"{line.strip()[:80]}"
                )
                continue
            row = AnnotationRow(
                seqname=seqname,
                source=fields[1].replace(" ", "_"),
                feature=fields[2],
                start=start,
                end=end,
                score=fields[5],
                strand=fields[6],
                frame=fields[7],
                attributes=parse_annot(fields[8]),
            )
            rows.append(row)
        self._all_rows = rows
        self._group()

    def _group(self) -> None:
        self.genes.clear()
        self.unassigned_rows.clear()
        for row in self._all_rows:
            gid = row.attributes.get("gene_id")
            if gid is None:
                self.unassigned_rows.append(row)
                continue
            gene = self.genes.get(gid)
            if gene is None:
                gene = Gene(gene_id=gid)
                self.genes[gid] = gene
            gene.rows.append(row)
            if row.feature == "gene":
                continue
            tid = row.attributes.get("transcript_id")
            if tid is None:
                continue
            tx = gene.transcripts.get(tid)
            if tx is None:
                tx = Transcript(gene_id=gid, transcript_id=tid)
                gene.transcripts[tid] = tx
            tx.rows.append(row)

    def prune(
        self,
        *,
        valid_coords: bool = True,
        keep_transcript_ctx: bool = True,
    ) -> int:
        """Drop rows that a GTF-cleaning workflow cannot use.

        Applies the shared policy: rows with invalid coordinates, rows without
        ``gene_id``, and non-gene rows without ``transcript_id`` are removed
        (count returned). A ``transcript_id`` accidentally present on a ``gene``
        row is stripped. ``iter_rows()`` then only sees surviving rows.
        """
        survivors = []
        removed = 0
        for row in self._all_rows:
            if valid_coords and (row.start < 1 or row.end < row.start):
                LOGGER.warning(f"Skipping row with invalid coordinates: {row}")
                removed += 1
                continue
            if row.attributes.get("gene_id") is None:
                LOGGER.warning(f"Skipping row without gene_id: {row}")
                removed += 1
                continue
            if (
                keep_transcript_ctx
                and row.feature != "gene"
                and "transcript_id" not in row.attributes
            ):
                LOGGER.warning(f"Skipping non-gene row without transcript_id: {row}")
                removed += 1
                continue
            if row.feature == "gene":
                row.attributes.pop("transcript_id", None)
            survivors.append(row)
        self._all_rows = survivors
        self._group()
        return removed

    # -- access ------------------------------------------------------------

    def gene(self, gene_id: str) -> Optional[Gene]:
        return self.genes.get(gene_id)

    def iter_genes(self) -> Iterator[Gene]:
        yield from self.genes.values()

    def iter_transcripts(self) -> Iterator[Transcript]:
        for gene in self.genes.values():
            yield from gene.transcripts.values()

    def iter_rows(self) -> Iterator[AnnotationRow]:
        """Every parsed row in file order (incl. unassigned ones)."""
        yield from self._all_rows

    def remove_gene(self, gene_id: str) -> None:
        """Drop a gene (and its rows) from the annotation."""
        self.genes.pop(gene_id, None)
        self._all_rows = [
            r for r in self._all_rows if r.attributes.get("gene_id", None) != gene_id
        ]

    # -- output -------------------------------------------------------------

    def write_gtf(
        self,
        path: str,
        *,
        sort: bool = True,
        check: bool = False,
        extra_comments: tuple[str, ...] = (),
        bgzip: bool = True,
    ) -> None:
        """Serialize the annotation back to a GTF file.

        Attributes are written in insertion order; when ``sort`` is set the
        body is sorted by (seqname, start, end). ``bgzip`` additionally creates
        ``path + ".gz"`` with a tabix index (best-effort, using pysam).
        """
        with open(path, "w") as fw:
            for comment in self.header:
                fw.write(comment + "\n")
            for comment in extra_comments:
                comment = comment.strip()
                fw.write(
                    ("#" + comment if not comment.startswith("#") else comment) + "\n"
                )
            body = []
            for gene in self.genes.values():
                for row in gene.rows:
                    if check:
                        self._check_row(row)
                    body.append(row.to_gtf_line())
            if sort:
                body.sort(
                    key=lambda line: (
                        line.split("\t")[0],
                        int(line.split("\t")[3]),
                        int(line.split("\t")[4]),
                    )
                )
            fw.writelines(line + "\n" for line in body)

        if bgzip:
            try:
                gz = path + ".gz"
                pysam.tabix_compress(path, gz, force=True)
                pysam.tabix_index(gz, preset="gff", force=True)
            except Exception as e:  # pragma: no cover - best-effort indexing
                LOGGER.warning(f"Skipping bgzip/tabix indexing: {e}")

    @staticmethod
    def _check_row(row: AnnotationRow) -> None:
        """Invariants the refined output guarantees (strict mode)."""
        if "gene_id" not in row.attributes:
            raise ValueError(f"Refined row has no gene_id: {row}")
        if "gene_name" not in row.attributes:
            raise ValueError(f"Refined row has no gene_name: {row}")
        if "gene_type" not in row.attributes:
            raise ValueError(f"Refined row has no gene_type: {row}")
        if row.feature == "gene":
            if "transcript_id" in row.attributes:
                raise ValueError(f"Gene row has a transcript_id: {row}")
        else:
            for key in ("transcript_id", "transcript_name", "transcript_type"):
                if key not in row.attributes:
                    raise ValueError(f"Refined row has no {key}: {row}")
