#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Tests for the migrated metagene subcommand (`coralsnake metagene`).
# Ported from the standalone `metagene` package into coralsnake.
#
# These use the R64-1-1 yeast GTF bundled in tests/data and build an input
# sites file from its exons so the tests are self-contained (no network).

from pathlib import Path

import pytest

pytest.importorskip(
    "coralsnake.gtf", reason="metagene deps not available", exc_type=ImportError
)

from coralsnake.annotation import map_to_transcripts, normalize_positions  # noqa: E402
from coralsnake.gtf import load_gtf, prepare_exon_ref  # noqa: E402
from coralsnake.io import load_sites  # noqa: E402
from coralsnake.map_to_local import map_to_local  # noqa: E402

DATA = Path(__file__).resolve().parent / "data"
GTF = DATA / "R64-1-1.release57.gtf"


def _build_sites(tmp_path, n=120):
    """Build a small sites TSV from exon midpoints on chromosome I."""
    sites = []
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, feature = parts[0], parts[2]
            if chrom == "I" and feature == "exon":
                start, end = int(parts[3]), int(parts[4])
                sites.append(f"I\t{(start + end) // 2}\t{parts[6]}")
                if len(sites) >= n:
                    break
    path = tmp_path / "sites.tsv"
    path.write_text("ref\tposition\tstrand\n" + "\n".join(sites) + "\n")
    return str(path)


class TestGtf:
    def test_prepare_exon_ref_schema(self):
        ref = prepare_exon_ref(str(GTF))
        for col in [
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "transcript_id",
            "gene_id",
            "Start_exon",
            "End_exon",
            "start_codon_pos",
            "stop_codon_pos",
            "transcript_length",
            "exon_number",
        ]:
            assert col in ref.columns, f"missing column {col}"

    def test_load_gtf_caches(self):
        ref = load_gtf(str(GTF))
        assert ref.height > 0


class TestMetagenePipeline:
    def test_map_to_transcripts(self, tmp_path):
        ref = load_gtf(str(GTF))
        sites = load_sites(
            _build_sites(tmp_path), with_header=True, meta_col_index=[0, 1, 2]
        )
        annotated = map_to_transcripts(sites, ref)
        # All sites are on exons so every site should be annotated
        assert annotated.height == sites.height
        assert annotated["gene_id"].is_not_null().all()
        # 0-based half-open coordinates
        assert (annotated["Start"] >= 0).all()

    def test_normalize_positions(self, tmp_path):
        ref = load_gtf(str(GTF))
        sites = load_sites(
            _build_sites(tmp_path), with_header=True, meta_col_index=[0, 1, 2]
        )
        annotated = map_to_transcripts(sites, ref)
        gene_bins, gene_stats, gene_splits = normalize_positions(
            annotated, split_strategy="median", bin_number=100
        )
        assert gene_bins.height == 100
        assert len(gene_splits) == 3
        # Splits sum to ~1
        assert round(sum(gene_splits), 6) == 1.0
        assert "count" in gene_bins.columns


class TestMapToLocal:
    """Sanity tests for the migrated map_to_local with real exon data."""

    def test_real_exon_data(self, tmp_path):
        ref = load_gtf(str(GTF))
        sites = load_sites(
            _build_sites(tmp_path), with_header=True, meta_col_index=[0, 1, 2]
        )
        result = map_to_local(
            sites, ref, ref_id_col="transcript_id", keep_global_loc=True
        )
        assert result.height > 0
        assert "Chromosome" in result.columns
        assert (result["End"] > result["Start"]).all()
        assert "Start_global" in result.columns
