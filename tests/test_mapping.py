"""Tests for coralsnake.mapping – integration mapping."""

from pathlib import Path

import pytest
import pysam

mapping = pytest.importorskip(
    "coralsnake.mapping", reason="mapping deps (bwamem) not available"
)

score_to_mapq = mapping.score_to_mapq


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
            output_files=[out],
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
            output_files=[out],
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
            output_files=[out],
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

    def test_pe_flags_correctness(self, tmp_path):
        """Verify that PE reads have correct flags when one is RCed."""
        ref_fa = tmp_path / "repro_ref.fa"
        ref_seq = "CGCTCCTTCGGTGCTCTTGGCTGGGTGTCCCGCGGGGCCCGGGGCGT"
        with open(ref_fa, "w") as f:
            f.write(f">ref\n{ref_seq}{'A' * 100}\n")
        
        r1_fq = tmp_path / "repro_R1.fq"
        with open(r1_fq, "w") as f:
            f.write(f"@read1\n{ref_seq}\n+\n{'I' * len(ref_seq)}\n")
        
        # R2 is the RC of R1
        from coralsnake.utils import reverse_complement
        r2_seq = reverse_complement(ref_seq)
        r2_fq = tmp_path / "repro_R2.fq"
        with open(r2_fq, "w") as f:
            f.write(f"@read1\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")
        
        out_bam = tmp_path / "repro_out.bam"
        idx_dir = tmp_path / "repro_idx"
        
        mapping.map_file(
            r1_file=str(r1_fq),
            r2_file=str(r2_fq),
            ref_files=[str(ref_fa)],
            output_files=[str(out_bam)],
            index_dir=str(idx_dir),
            forward_library=True,
            min_alignment_length=10,
            min_mapping_ratio=0.8,
            threads=1
        )
        
        with pysam.AlignmentFile(str(out_bam), "rb") as bam:
            reads = list(bam)
            assert len(reads) == 2
            r1 = [r for r in reads if r.is_read1][0]
            r2 = [r for r in reads if r.is_read2][0]
            
            # R1 should be forward (Flag 99 = 1+2+32+64)
            # 0x1: paired, 0x2: proper pair, 0x20: mate reverse, 0x40: read1
            assert not r1.is_reverse
            assert r1.mate_is_reverse
            assert r1.flag == 99
            
            # R2 should be reverse (Flag 147 = 1+2+16+128)
            # 0x1: paired, 0x2: proper pair, 0x10: reverse, 0x80: read2
            assert r2.is_reverse
            assert not r2.mate_is_reverse
            assert r2.flag == 147

    def test_index_only(self, tmp_path, data_dir):
        if not (data_dir / "ref.fa").exists():
            pytest.skip("ref.fa not present")
        idx = tmp_path / "idx"
        mapping.map_file(
            r1_file=None,
            r2_file=None,
            ref_files=[str(data_dir / "ref.fa")],
            output_files=[str(tmp_path / "dummy.bam")],
            index_dir=str(idx),
            index_only=True,
        )
        assert (idx / "ref.orig.fa").exists()
        assert (idx / "ref.mk.amb").exists()
