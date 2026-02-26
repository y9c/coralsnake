"""Tests for coralsnake.mapping – scoring helpers and integration mapping."""

from pathlib import Path

import pytest
import pysam

mapping = pytest.importorskip("coralsnake.mapping", reason="mapping deps (bwamem) not available")

score_to_mapq = mapping.score_to_mapq
cal_md_and_tag = mapping.cal_md_and_tag
calculate_directional_score = mapping.calculate_directional_score
find_properly_paired_hits = mapping.find_properly_paired_hits
filter_hits = mapping.filter_hits


# ---------------------------------------------------------------------------
# score_to_mapq
# ---------------------------------------------------------------------------
class TestScoreToMapq:
    def test_clamped_high(self):
        assert score_to_mapq(100) == 60

    def test_clamped_low(self):
        assert score_to_mapq(-5) == 0

    def test_passthrough(self):
        assert score_to_mapq(30) == 30


# ---------------------------------------------------------------------------
# cal_md_and_tag (delegates to C)
# ---------------------------------------------------------------------------
class TestCalMdAndTag:
    def test_delegates(self):
        md, *rest = cal_md_and_tag("4M", "ACGT", "ACGT", True)
        assert md == "4"


# ---------------------------------------------------------------------------
# calculate_directional_score (delegates to C)
# ---------------------------------------------------------------------------
class TestCalculateDirectionalScore:
    def test_delegates(self):
        score, w, b = calculate_directional_score("4M", "ACGT", "ACGT", True)
        assert score == 4


# ---------------------------------------------------------------------------
# find_properly_paired_hits
# ---------------------------------------------------------------------------
class TestFindProperlyPairedHits:
    class FakeHit:
        def __init__(self, ctg, r_st, r_en, strand, read_num):
            self.ctg = ctg
            self.r_st = r_st
            self.r_en = r_en
            self.strand = strand
            self.read_num = read_num

    def test_close_pair(self):
        h1 = self.FakeHit("chr1", 100, 200, 1, 1)
        h2 = self.FakeHit("chr1", 300, 400, 1, 2)
        pairs = find_properly_paired_hits([h1, h2])
        assert len(pairs) == 1
        assert pairs[0] == (h1, h2)

    def test_too_far(self):
        h1 = self.FakeHit("chr1", 0, 100, 1, 1)
        h2 = self.FakeHit("chr1", 5000, 5100, 1, 2)
        pairs = find_properly_paired_hits([h1, h2])
        assert len(pairs) == 0

    def test_different_contigs(self):
        h1 = self.FakeHit("chr1", 100, 200, 1, 1)
        h2 = self.FakeHit("chr2", 100, 200, 1, 2)
        pairs = find_properly_paired_hits([h1, h2])
        assert len(pairs) == 0

    def test_same_read_num_not_paired(self):
        h1 = self.FakeHit("chr1", 100, 200, 1, 1)
        h2 = self.FakeHit("chr1", 300, 400, 1, 1)
        pairs = find_properly_paired_hits([h1, h2])
        assert len(pairs) == 0

    def test_symmetric_distance(self):
        h1 = self.FakeHit("chr1", 500, 600, 1, 1)
        h2 = self.FakeHit("chr1", 100, 200, 1, 2)
        pairs = find_properly_paired_hits([h1, h2])
        assert len(pairs) == 1


# ---------------------------------------------------------------------------
# filter_hits
# ---------------------------------------------------------------------------
class TestFilterHits:
    class FakeHit:
        def __init__(self, mapq, blen, mlen, read_num):
            self.mapq = mapq
            self.blen = blen
            self.mlen = mlen
            self.read_num = read_num

    def test_passes(self):
        h = self.FakeHit(mapq=10, blen=50, mlen=40, read_num=1)
        result = filter_hits([h], "A" * 50, None, min_alignment_length=20, min_mapping_ratio=0.5)
        assert len(result) == 1

    def test_low_mapq(self):
        h = self.FakeHit(mapq=0, blen=50, mlen=40, read_num=1)
        result = filter_hits([h], "A" * 50, None)
        assert len(result) == 0

    def test_short_alignment(self):
        h = self.FakeHit(mapq=10, blen=10, mlen=10, read_num=1)
        result = filter_hits([h], "A" * 50, None, min_alignment_length=20)
        assert len(result) == 0

    def test_low_ratio(self):
        h = self.FakeHit(mapq=10, blen=50, mlen=10, read_num=1)
        result = filter_hits([h], "A" * 50, None, min_mapping_ratio=0.5)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# Integration: map_file (requires bwamem at runtime)
# ---------------------------------------------------------------------------
class TestMapFileIntegration:
    @pytest.fixture(autouse=True)
    def _check_bwamem(self):
        pytest.importorskip("bwamem", reason="bwamem not installed")

    def test_pe_forward(self, tmp_path, data_dir, has_integration_data):
        out = str(tmp_path / "out.bam")
        idx = str(tmp_path / "idx")
        mapping.map_file(
            r1_file=str(data_dir / "test1.fq"),
            r2_file=str(data_dir / "test2.fq"),
            ref_files=[str(data_dir / "ref.fa")],
            output_files=out,
            forward_library=True,
            max_mismatches=0,
            threads=2,
            min_alignment_length=8,
            min_mapping_ratio=0.5,
            index_dir=idx,
            batch_size=50,
        )
        with pysam.AlignmentFile(out, "rb") as f:
            assert sum(1 for _ in f) > 0

    def test_pe_reverse(self, tmp_path, data_dir, has_integration_data):
        out = str(tmp_path / "out.bam")
        idx = str(tmp_path / "idx")
        mapping.map_file(
            r1_file=str(data_dir / "test1.fq"),
            r2_file=str(data_dir / "test2.fq"),
            ref_files=[str(data_dir / "ref.fa")],
            output_files=out,
            forward_library=False,
            max_mismatches=0,
            threads=2,
            min_alignment_length=8,
            min_mapping_ratio=0.5,
            index_dir=idx,
            batch_size=50,
        )
        with pysam.AlignmentFile(out, "rb") as f:
            assert sum(1 for _ in f) > 0

    def test_se_with_unmap(self, tmp_path, data_dir, has_integration_data):
        out = str(tmp_path / "out.bam")
        unmap = str(tmp_path / "unmap.bam")
        idx = str(tmp_path / "idx")
        mapping.map_file(
            r1_file=str(data_dir / "test1.fq"),
            r2_file=None,
            ref_files=[str(data_dir / "ref.fa")],
            output_files=out,
            unmap_file=unmap,
            forward_library=True,
            max_mismatches=0,
            threads=2,
            min_alignment_length=8,
            min_mapping_ratio=0.5,
            index_dir=idx,
            batch_size=50,
        )
        assert Path(out).exists()
        assert Path(unmap).exists()

    def test_index_only(self, tmp_path, data_dir):
        if not (data_dir / "ref.fa").exists():
            pytest.skip("ref.fa not present")
        idx = tmp_path / "idx"
        mapping.map_file(
            r1_file=None,
            r2_file=None,
            ref_files=[str(data_dir / "ref.fa")],
            output_files=str(tmp_path / "dummy.bam"),
            index_dir=str(idx),
            index_only=True,
        )
        assert (idx / "ref.orig.fa").exists()
        assert (idx / "ref.mk.amb").exists()
