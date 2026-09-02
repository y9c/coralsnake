#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Convert a built-in (or custom) exon reference frame - the
# ``prepare_exon_ref`` schema, as served by the ``data`` release and cached
# under ``~/.cache/coralsnake/`` - into the other input formats used across
# the pipeline, so a single small download serves every tool:
#
#   * the ``prepare`` annotation table (TSV) -> ``liftover -a`` /
#     ``liftover --table`` / ``annotate --annotation``
#   * a GTF -> ``annotate --reference-gtf`` (full variant effect)
#
# Coordinate conventions: the reference stores 0-based half-open [Start, End)
# genomic positions and 5'->3' spliced offsets (Start_exon/End_exon). The
# exports convert back to the 1-based inclusive conventions each consumer
# expects (``Transcript.to_tsv`` / standard GTF).

import polars as pl

from .utils import get_logger

LOGGER = get_logger(__name__)

# Columns every reference must have to be usable (v1 core schema; v2 adds
# gene_name / gene_biotype, which are optional here).
_REQUIRED = (
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "gene_id",
    "transcript_id",
    "transcript_length",
    "Start_exon",
    "End_exon",
)


def _check_schema(df: pl.DataFrame) -> None:
    missing = [c for c in _REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(
            f"reference is missing required columns {missing}; the cached file "
            "may be corrupt or from an incompatible version - delete it and "
            "re-download (coralsnake reference download <ref>)"
        )


def _codon_1based(pcol: str) -> pl.Expr:
    """1-based genomic coordinate of the first base (5' end) of a codon.

    ``pcol`` is the codon's 0-based position in the spliced transcript.
    Evaluated per exon row: only the row whose exon actually contains the
    codon (``Start_exon <= p < End_exon``) yields a value; every other row of
    the transcript yields null, so ``first()`` over the transcript group
    recovers the position exactly once.
    """
    p = pl.col(pcol)
    return (
        pl.when(
            p.is_not_null() & (p >= pl.col("Start_exon")) & (p < pl.col("End_exon"))
        )
        .then(
            pl.when(pl.col("Strand") == "+")
            .then(pl.col("Start") + (p - pl.col("Start_exon")))
            .otherwise(pl.col("End") - 1 - (p - pl.col("Start_exon")))
            + 1
        )
        .otherwise(None)
    )


def export_table(df: pl.DataFrame, path: str) -> str:
    """Write the ``prepare``-format annotation table (TSV) for ``df``.

    The layout matches ``coralsnake prepare -c -n -t -x`` (same header order)
    so the file can be fed directly to ``liftover -a``, ``liftover --table``
    and ``annotate --annotation``. ``spans`` are 1-based inclusive in 5'->3'
    transcript order; the codon columns hold the 1-based genomic coordinate
    of the codon's first base (empty when the transcript has no codon).
    """
    _check_schema(df)
    has_name = "gene_name" in df.columns
    has_bio = "gene_biotype" in df.columns

    # span strings: "start-end" (1-based inclusive), one per exon row
    work = df.with_columns(
        _codon_1based("start_codon_pos").alias("_sc"),
        _codon_1based("stop_codon_pos").alias("_sp"),
        ((pl.col("Start") + 1).cast(pl.Utf8) + "-" + pl.col("End").cast(pl.Utf8)).alias(
            "_span"
        ),
    )

    aggs = [
        pl.col("Chromosome").first().alias("chrom"),
        pl.col("Strand").first().alias("strand"),
        pl.col("_span").str.join(",").alias("spans"),
        # max() is null-ignoring: exactly one row per transcript carries the
        # codon position (by construction of _codon_1based), so this recovers
        # it no matter which exon the codon lies in.
        pl.col("_sc").max().alias("start_codon"),
        pl.col("_sp").max().alias("stop_codon"),
        pl.col("Start").min().alias("_tx_start"),
        pl.col("End").max().alias("_tx_end"),
    ]
    if has_name:
        aggs.append(pl.col("gene_name").max().alias("gene_name"))
    if has_bio:
        aggs.append(pl.col("gene_biotype").max().alias("transcript_biotype"))

    grouped = (
        work.sort(["gene_id", "transcript_id", "Start_exon"])
        .group_by(["gene_id", "transcript_id"], maintain_order=True)
        .agg(aggs)
        .with_columns(
            (pl.col("_tx_start") + 1).cast(pl.Utf8).alias("transcript_start"),
            pl.col("_tx_end").cast(pl.Utf8).alias("transcript_end"),
        )
        .drop(["_tx_start", "_tx_end"])
    )

    cols = [
        "gene_id",
        "transcript_id",
        "chrom",
        "strand",
        "spans",
        "start_codon",
        "stop_codon",
    ]
    if has_name:
        cols.append("gene_name")
    if has_bio:
        cols.append("transcript_biotype")
    cols += ["transcript_start", "transcript_end"]

    # Codon cells that are null stay null: the CSV writer emits null as an
    # empty field, which readers treat as "no codon" (a literal empty string
    # would be quoted as '""' and read back as data).
    grouped.select(cols).write_csv(path, separator="\t", include_header=True)
    LOGGER.info(f"[green]✓ Wrote annotation table to {path}[/green]")
    return path


def export_gtf(df: pl.DataFrame, path: str) -> str:
    """Re-emit ``df`` as a GTF (gene / transcript / exon / codon rows).

    The output is parseable by ``prepare_exon_ref`` (it round-trips every
    stable column, including codon positions re-derived from the emitted
    genomic coordinates) and is a valid input for
    ``annotate --reference-gtf``. GTF conventions: 1-based inclusive,
    coordinates ascending regardless of strand, ``source`` = ``Ensembl``.
    Exon rows carry ``exon_number``; canonical transcripts
    (``transcript_level == 0``) carry ``tag "Ensembl_canonical"``; other
    levels carry ``transcript_support_level`` (6/10 - untrusted - are
    omitted, mirroring how the level was derived).
    """
    _check_schema(df)
    has_name = "gene_name" in df.columns
    has_bio = "gene_biotype" in df.columns

    work = df.with_columns(
        _codon_1based("start_codon_pos").alias("_sc1"),
        _codon_1based("stop_codon_pos").alias("_sp1"),
    )

    tx = (
        work.group_by(["gene_id", "transcript_id"], maintain_order=True)
        .agg(
            pl.col("Chromosome").first(),
            pl.col("Strand").first(),
            pl.col("Start").min().alias("tx_start"),
            pl.col("End").max().alias("tx_end"),
            pl.col("transcript_level").max(),
            # max() is null-ignoring (exactly one row per transcript is
            # non-null by construction of _codon_1based).
            pl.col("_sc1").max().alias("start_codon_g"),
            pl.col("_sp1").max().alias("stop_codon_g"),
            *([pl.col("gene_name").max()] if has_name else []),
            *([pl.col("gene_biotype").max()] if has_bio else []),
        )
        .sort(["gene_id", "transcript_id"])
        .to_dicts()
    )

    gene_bbox = {
        r["gene_id"]: (int(r["g_start"]) + 1, int(r["g_end"]))
        for r in work.group_by("gene_id")
        .agg(
            pl.col("Start").min().alias("g_start"),
            pl.col("End").max().alias("g_end"),
        )
        .to_dicts()
    }

    # exon rows per transcript, in 5'->3' order (Start_exon ascending)
    exons_by_tx: dict[str, list[dict]] = {}
    for r in (
        work.sort(["gene_id", "transcript_id", "Start_exon"])
        .select(
            [
                "transcript_id",
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "exon_number",
            ]
        )
        .to_dicts()
    ):
        exons_by_tx.setdefault(r["transcript_id"], []).append(r)

    def _attrs(mapping: dict) -> str:
        return (
            " ".join(f'{k} "{v}"' for k, v in mapping.items() if v not in (None, ""))
            + ";"
        )

    lines: list[str] = []
    seen_genes: set[str] = set()
    for t in tx:
        chrom, strand = t["Chromosome"], t["Strand"]
        ts, te = int(t["tx_start"]) + 1, int(t["tx_end"])
        level = t["transcript_level"]

        if t["gene_id"] not in seen_genes:
            seen_genes.add(t["gene_id"])
            gs, ge = gene_bbox[t["gene_id"]]
            ga = {"gene_id": t["gene_id"]}
            if has_name and t.get("gene_name"):
                ga["gene_name"] = t["gene_name"]
            if has_bio and t.get("gene_biotype"):
                ga["gene_biotype"] = t["gene_biotype"]
            lines.append(
                f"{chrom}\tEnsembl\tgene\t{gs}\t{ge}\t.\t{strand}\t.\t{_attrs(ga)}"
            )

        ta = {"gene_id": t["gene_id"]}
        if has_name and t.get("gene_name"):
            ta["gene_name"] = t["gene_name"]
        if has_bio and t.get("gene_biotype"):
            ta["gene_biotype"] = t["gene_biotype"]
        ta["transcript_id"] = t["transcript_id"]
        exons = exons_by_tx.get(t["transcript_id"], [])
        ta["exon_count"] = str(len(exons))
        if level == 0:
            ta["tag"] = "Ensembl_canonical"
        elif isinstance(level, (int, float)) and level not in (6, 10):
            ta["transcript_support_level"] = str(int(level))
        lines.append(
            f"{chrom}\tEnsembl\ttranscript\t{ts}\t{te}\t.\t{strand}\t.\t{_attrs(ta)}"
        )

        for i, e in enumerate(sorted(exons, key=lambda e: e["Start"])):
            es, ee = int(e["Start"]) + 1, int(e["End"])
            ea = {
                "gene_id": t["gene_id"],
                "transcript_id": t["transcript_id"],
                "exon_number": e["exon_number"] or str(i + 1),
            }
            # gene_name/gene_biotype on exon rows too (modern Ensembl GTFs
            # carry them on every row) so a re-parse round-trips the v2
            # identity columns.
            if has_name and t.get("gene_name"):
                ea["gene_name"] = t["gene_name"]
            if has_bio and t.get("gene_biotype"):
                ea["gene_biotype"] = t["gene_biotype"]
            lines.append(
                f"{e['Chromosome']}\tEnsembl\texon\t{es}\t{ee}\t.\t{e['Strand']}\t.\t{_attrs(ea)}"
            )

        for feature, g1 in (
            ("start_codon", t["start_codon_g"]),
            ("stop_codon", t["stop_codon_g"]),
        ):
            if g1 is None:
                continue
            g1 = int(g1)
            # 1-based inclusive codon bounds in GTF (ascending) convention:
            # on the minus strand the 5' (first) base is the HIGHEST genomic
            # coordinate, so the codon occupies [g1-2, g1] inclusive.
            cs, ce = (g1, g1 + 2) if strand == "+" else (g1 - 2, g1)
            ca = {"gene_id": t["gene_id"], "transcript_id": t["transcript_id"]}
            lines.append(
                f"{chrom}\tEnsembl\t{feature}\t{cs}\t{ce}\t.\t{strand}\t.\t{_attrs(ca)}"
            )

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    LOGGER.info(f"[green]✓ Wrote GTF to {path}[/green]")
    return path
