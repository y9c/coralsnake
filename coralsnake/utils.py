#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-25 14:21

import logging
import os
import textwrap
from collections import defaultdict
from pathlib import Path

import pysam
from rich.console import Console
from rich.logging import RichHandler

from . import seqops

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


def require_plotting(backend: str | None = "Agg"):
    """Import matplotlib lazily, raising a helpful error if it is missing.

    matplotlib is an optional dependency (heavy); users should install it with
    ``pip install 'coralsnake[plot]'``. ``backend`` (default ``"Agg"``) is
    applied before pyplot is imported; pass ``None`` to leave the backend alone.
    """
    try:
        import matplotlib
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "Plotting requires the optional 'matplotlib' dependency.\n"
            "Install it with:  pip install 'coralsnake[plot]'"
        ) from e
    if backend is not None:
        matplotlib.use(backend)
    return matplotlib


def reverse_complement(seq: str) -> str:
    """Fast reverse complement using C implementation."""
    return seqops.reverse_complement(seq)


def interval_groups(labels, labels2=None):
    """Map categorical string labels to shared ``uint32`` group ids for ruranges.

    ruranges' ``overlaps``/``groups=``/``groups2=`` arguments require a
    ``uint32`` id per row, with equal labels sharing an id. This centralises that
    (previously copy-pasted) conversion: labels and (optionally) labels2 are
    numbered together so identical labels map to the same id. Fully vectorized
    via ``np.unique(..., return_inverse=True)`` (same deterministic ordering the
    callers previously used, but with vectorized conversion).
    """
    import numpy as np

    a = np.asarray(labels, dtype=str)
    if labels2 is None:
        _, inv = np.unique(a, return_inverse=True)
        return inv.astype(np.uint32)

    b = np.asarray(labels2, dtype=str)
    both = np.concatenate([a, b])
    _, inv = np.unique(both, return_inverse=True)
    return inv[: a.size].astype(np.uint32), inv[a.size :].astype(np.uint32)


def convert_file_realtime(
    input_file: str,
    output_file: str,
    base_from: str,
    base_to: str,
):
    """
    Convert FASTA file using fast C implementation (line-by-line, memory efficient).

    This function is optimized for large genome reference files. It reads and writes
    line-by-line without loading the entire file into memory, making it much faster
    and more memory-efficient than dnaio-based conversion.

    Args:
        input_file: Path to input FASTA file
        output_file: Path to output converted FASTA file
        base_from: Bases to convert from (e.g., "AC" for MK conversion)
        base_to: Bases to convert to (e.g., "GT" for MK conversion)

    Examples:
        MK conversion (A->G, C->T): convert_file_realtime(in, out, "AC", "GT")
        KM conversion (G->A, T->C): convert_file_realtime(in, out, "GT", "AC")

    Note:
        - Only works with FASTA format (not FASTQ)
        - Does not add YS:Z tags
        - Header lines (starting with '>') are written directly
        - Sequence lines are converted base-by-base
    """
    seqops.convert_fasta_file(input_file, output_file, base_from, base_to)


class Span:
    def __init__(self, start: int, end: int):
        # 0-based
        self.start = start
        self.end = end

    def __iter__(self):
        yield self.start
        yield self.end

    def __repr__(self) -> str:
        return f"(start={self.start}, end={self.end})"


class Transcript:
    def __init__(
        self,
        gene_id: str | None = None,
        transcript_id: str | None = None,
        chrom: str = "",
        strand: str = "",
        exons: dict[str | int, Span] | None = None,
        gene_name: str | None = None,
    ):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.exons = {} if exons is None else exons
        # calculated feature
        self._exons_forwards = None
        self._cum_exon_lens = None
        self._seq = None
        # extra feature
        self.gene_name: str | None = gene_name
        self.transcript_biotype: str | None = None
        self.start_codon: Span | None = None
        self.stop_codon: Span | None = None
        self.priority: tuple[int, int] = (10, 0)
        self.conflict: bool = False
        self._annotations: dict = {}

    def add_exon(self, exon_id: str | int, exon: Span) -> None:
        self.exons[exon_id] = exon
        self._exons_forwards = None
        self._cum_exon_lens = None
        self._seq = None

    def sort_exons(self) -> None:
        self.exons = dict(
            sorted(
                self.exons.items(), key=lambda x: x[1].start, reverse=self.strand == "-"
            ),
        )
        self._exons_forwards = None
        self._cum_exon_lens = None
        self._seq = None

    def to_tsv(
        self,
        with_codon=False,
        with_genename=False,
        with_biotype=False,
        with_txpos=False,
    ) -> str:
        # convert into 1-based
        line = []
        for key in ["gene_id", "transcript_id", "chrom", "strand"]:
            v = getattr(self, key)
            if v is not None:
                line.append(v)
            else:
                line.append("")
        # spans are emitted in 5'->3' transcript order (deterministic,
        # independent of the exons dict's insertion order, which get_seq()
        # may re-sort) - table-mode annotation accumulates offsets in this
        # order
        exons = sorted(self.exons.values(), key=lambda e: e.start)
        if self.strand == "-":
            exons.reverse()
        line += [",".join([f"{v.start + 1}-{v.end}" for v in exons])]
        if with_codon:
            # use the start postion of the span
            line += [
                str(self.start_codon.start + 1) if self.start_codon else "",
                str(self.stop_codon.start + 1) if self.stop_codon else "",
            ]
        if with_genename:
            line.append(self.gene_name if self.gene_name is not None else "")
        if with_biotype:
            line.append(self.transcript_biotype if self.transcript_biotype else "")
        # transcript start and end: the genomic bounding box (1-based)
        # [min_exon_start, max_exon_end], identical for both strands.
        if with_txpos:
            es = list(self.exons.values())
            line += [
                str(min(e.start for e in es) + 1),
                str(max(e.end for e in es)),
            ]

        return "\t".join(line)

    def get_genome_spans(self) -> str:
        # 1-based
        return ",".join([f"{v.start + 1}-{v.end}" for v in self.exons.values()])

    def get_gene_spans(self) -> str:
        # 1-based
        gene_spans = []
        exon_start = 1
        for exon_end in self.cum_exon_lens:
            gene_spans.append(f"{exon_start}-{exon_end}")
            exon_start = exon_end + 1
        return ",".join(gene_spans)

    def get_seq(self, fasta: pysam.FastaFile, sort=True, upper=True, wrap=0):
        if self._seq is None:
            if sort:
                self.sort_exons()
            seq = ""
            for _, v in self.exons.items():
                e = fasta.fetch(self.chrom, v.start, v.end)
                if self.strand == "-":
                    e = reverse_complement(e)
                seq += e
            self._seq = seq

        seq_out = self._seq
        if upper:
            seq_out = seq_out.upper()
        if wrap > 0:
            seq_out = textwrap.fill(seq_out, wrap)
        return seq_out

    @property
    def exons_forwards(self) -> list[Span]:
        if self._exons_forwards is None:
            self._exons_forwards = sorted(self.exons.values(), key=lambda e: e.start)
        return self._exons_forwards

    @property
    def cum_exon_lens(self) -> list[int]:
        if self._cum_exon_lens is None:
            lengths = [exon.end - exon.start for exon in self.exons_forwards]
            cum_lengths = []
            total = 0
            for length in lengths:
                total += length
                cum_lengths.append(total)
            self._cum_exon_lens = cum_lengths
        return self._cum_exon_lens

    @property
    def length(self) -> int:
        return (self.cum_exon_lens or [0])[-1]

    def __repr__(self) -> str:
        res = []
        for key in [
            "gene_id",
            "transcript_id",
            "chrom",
            "strand",
            # "exons_forwards",
        ]:
            res.append(f"{key}={getattr(self, key)}")
        return f"Transcript({', '.join(res)})"


def load_annotation(
    annotation_file: str, with_header: bool = True
) -> dict[str, dict[str, Transcript]]:
    annot = defaultdict(dict)
    with open(annotation_file, "r") as f:
        if with_header:
            header = f.readline().strip("\n").split("\t")
            gene_id_idx = header.index("gene_id")
            transcript_id_idx = header.index("transcript_id") if "transcript_id" in header else gene_id_idx
            chrom_idx = header.index("chrom")
            strand_idx = header.index("strand")
            spans_idx = header.index("spans")
            if "gene_name" in header:
                gene_name_idx = header.index("gene_name")
            else:
                gene_name_idx = None
        else:
            (
                gene_id_idx,
                transcript_id_idx,
                chrom_idx,
                strand_idx,
                spans_idx,
                gene_name_idx,
            ) = 0, 1, 2, 3, 4, None
        for line in f:
            fields = line.strip("\n").split("\t")
            if len(fields) < 5:
                continue
            gene_id = fields[gene_id_idx]
            transcript_id = fields[transcript_id_idx]
            chrom = fields[chrom_idx]
            strand = fields[strand_idx]
            spans = fields[spans_idx]
            gene_name = fields[gene_name_idx] if gene_name_idx is not None else None
            exons = {
                idx: Span(int(start) - 1, int(end))
                for idx, span in enumerate(spans.split(","), 1)
                for start, end in [span.split("-")]
            }
            transcript = Transcript(
                gene_id=gene_id,
                transcript_id=transcript_id,
                chrom=chrom,
                strand=strand,
                exons=exons,
                gene_name=gene_name,
            )
            annot[gene_id][transcript_id] = transcript
    return annot


def load_faidx(faidx_file: str) -> dict[str, int]:
    faidx = {}
    with open(faidx_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            faidx[fields[0]] = int(fields[1])
    return faidx


def format_duration(sec: float) -> str:
    """
    Format duration in seconds to human-readable string.

    Args:
        sec: Duration in seconds

    Returns:
        Formatted duration string (e.g., "12.34s", "5m 23s", "2h 15m 30s")
    """
    if sec < 60:
        return f"{sec:.2f}s"
    if sec < 3600:
        minutes = int(sec // 60)
        seconds = int(round(sec % 60))
        return f"{minutes}m {seconds}s"
    hours = int(sec // 3600)
    rem = sec % 3600
    minutes = int(rem // 60)
    seconds = int(round(rem % 60))
    return f"{hours}h {minutes}m {seconds}s"


class NewlineRichHandler(RichHandler):
    """Rich logging handler that adds newlines for better visual separation."""

    def __init__(
        self,
        console=None,
        rich_tracebacks=True,
        markup=True,
        show_time=True,
        show_path=False,
        show_level=True,
        enable_link_path=False,
        **kwargs,
    ):
        if console is None:
            console = Console(stderr=True)
        super().__init__(
            console=console,
            rich_tracebacks=rich_tracebacks,
            markup=markup,
            show_time=show_time,
            show_path=show_path,
            show_level=show_level,
            enable_link_path=enable_link_path,
            **kwargs,
        )

    def emit(self, record):
        try:
            message = self.format(record)
            self.console.print("\n" + message)
        except Exception:
            self.handleError(record)


def setup_rich_logger(
    name: str = "coralsnake", level: int = logging.INFO, console: Console | None = None
) -> logging.Logger:
    """
    Set up and configure a rich logger for the application.

    Args:
        name: Name of the logger (default: "coralsnake")
        level: Logging level (default: INFO)
        console: Rich console instance to use (optional)

    Returns:
        Configured logger instance with rich formatting
    """
    logger = logging.getLogger(name)
    logger.handlers = []  # Remove any existing handlers
    logger.addHandler(NewlineRichHandler(console=console))
    logger.setLevel(level)
    logger.propagate = False  # Prevent propagation to root logger
    return logger


def setup_logger(name: str = "coralsnake", level: int = logging.INFO) -> logging.Logger:
    """
    Set up and configure a logger for the application.

    Args:
        name: Name of the logger (default: "coralsnake")
        level: Logging level (default: INFO)

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)

    # Only add handlers if they haven't been added already
    if not logger.handlers:
        logger.setLevel(level)

        # Create console handler with formatting
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    logger.propagate = False  # Prevent double logging
    return logger


def get_cache_dir(app_name: str = "coralsnake") -> Path:
    """
    Get the cache directory path following XDG Base Directory Specification.

    Args:
        app_name: Name of the application (default: "coralsnake")

    Returns:
        Path object pointing to the cache directory
    """
    if "XDG_CACHE_HOME" in os.environ:
        cache_home = Path(os.environ["XDG_CACHE_HOME"])
    else:
        cache_home = Path.home() / ".cache"

    cache_dir = cache_home / app_name
    cache_dir.mkdir(parents=True, exist_ok=True)

    return cache_dir


def ensure_dir(path: str | Path) -> Path:
    """
    Ensure a directory exists, creating it if necessary.

    Args:
        path: Directory path to ensure exists

    Returns:
        Path object for the directory
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_file_hash(file_path: str | Path) -> str:
    """
    Calculate MD5 hash based on file path and first chunk of content.
    This is much faster than hashing the entire file, especially for large GTF files.

    Args:
        file_path: Path to the file

    Returns:
        MD5 hash as a hexadecimal string
    """
    import hashlib

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    hash_md5 = hashlib.md5()

    # Hash the absolute file path
    hash_md5.update(str(file_path.resolve()).encode("utf-8"))

    # Hash the first chunk (4096 * 8 = 32KB) of the file
    chunk_size = 4096 * 8
    with open(file_path, "rb") as f:
        first_chunk = f.read(chunk_size)
        hash_md5.update(first_chunk)

    return hash_md5.hexdigest()


def get_file_size(file_path: str | Path) -> int:
    """
    Get the size of a file in bytes.

    Args:
        file_path: Path to the file

    Returns:
        File size in bytes
    """
    return Path(file_path).stat().st_size


def format_file_size(size_bytes: int) -> str:
    """
    Format file size in bytes to human-readable format.

    Args:
        size_bytes: Size in bytes

    Returns:
        Human-readable file size string
    """
    size = float(size_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} PB"
