from pathlib import Path

import pytest
import pysam

from coralsnake.mapping import map_file


def test_map_pe_forward_library(tmp_path: Path):
    """Test paired-end mapping with forward library."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref = data_dir / "ref.fa"
    r1 = data_dir / "test1.fq"
    r2 = data_dir / "test2.fq"
    if not (ref.exists() and r1.exists() and r2.exists()):
        pytest.skip("integration test data not present")

    out_bam = str(tmp_path / "out.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        r1_file=str(r1),
        r2_file=str(r2),
        ref_files=[str(ref)],
        output_files=out_bam,
        unmap_file=None,
        forward_library=True,
        max_mismatches=0,
        threads=2,
        min_alignment_length=8,
        min_mapping_ratio=0.5,
        index_dir=index_dir,
        index_only=False,
        batch_size=50,
    )

    bam = pysam.AlignmentFile(out_bam, "rb")
    assert sum(1 for _ in bam) > 0
    bam.close()


def test_map_pe_reverse_library(tmp_path: Path):
    """Test paired-end mapping with reverse library."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref = data_dir / "ref.fa"
    r1 = data_dir / "test1.fq"
    r2 = data_dir / "test2.fq"
    if not (ref.exists() and r1.exists() and r2.exists()):
        pytest.skip("integration test data not present")

    out_bam = str(tmp_path / "out.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        r1_file=str(r1),
        r2_file=str(r2),
        ref_files=[str(ref)],
        output_files=out_bam,
        unmap_file=None,
        forward_library=False,  # Reverse library
        max_mismatches=0,
        threads=2,
        min_alignment_length=8,
        min_mapping_ratio=0.5,
        index_dir=index_dir,
        index_only=False,
        batch_size=50,
    )

    bam = pysam.AlignmentFile(out_bam, "rb")
    reads = list(bam)
    bam.close()
    assert len(reads) > 0


def test_map_pe_integration(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref = data_dir / "ref.fa"
    r1 = data_dir / "test1.fq"
    r2 = data_dir / "test2.fq"
    if not (ref.exists() and r1.exists() and r2.exists()):
        pytest.skip("integration test data not present")

    out_bam = str(tmp_path / "out.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        r1_file=str(r1),
        r2_file=str(r2),
        ref_files=[str(ref)],
        output_files=out_bam,
        unmap_file=None,
        forward_library=True,
        max_mismatches=0,
        threads=2,
        min_alignment_length=8,
        min_mapping_ratio=0.5,
        index_dir=index_dir,
        index_only=False,
        batch_size=50,
    )

    bam = pysam.AlignmentFile(out_bam, "rb")
    assert sum(1 for _ in bam) > 0
    bam.close()


def test_index_only(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref = data_dir / "ref.fa"
    if not ref.exists():
        pytest.skip("integration test data not present")

    index_dir = tmp_path / "idx"
    map_file(
        r1_file=None,
        r2_file=None,
        ref_files=[str(ref)],
        output_files=str(tmp_path / "dummy.bam"),
        unmap_file=None,
        index_dir=str(index_dir),
        index_only=True,
    )

    # BWA index should exist (ref.orig.fa + ref.mk.fa*.{amb,ann,bwt,pac,sa})
    assert (index_dir / "ref.orig.fa").exists()
    assert (index_dir / "ref.mk.amb").exists()


def test_map_pe_with_unmap_file(tmp_path: Path):
    """Test mapping with unmapped reads written to separate file."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref = data_dir / "ref.fa"
    r1 = data_dir / "test1.fq"
    r2 = data_dir / "test2.fq"
    if not (ref.exists() and r1.exists() and r2.exists()):
        pytest.skip("integration test data not present")

    out_bam = str(tmp_path / "out.bam")
    unmap_bam = str(tmp_path / "unmap.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        r1_file=str(r1),
        r2_file=str(r2),
        ref_files=[str(ref)],
        output_files=out_bam,
        unmap_file=unmap_bam,
        forward_library=True,
        max_mismatches=0,
        threads=2,
        min_alignment_length=8,
        min_mapping_ratio=0.5,
        index_dir=index_dir,
        index_only=False,
        batch_size=50,
    )

    # Both files should exist
    assert Path(out_bam).exists()
    assert Path(unmap_bam).exists()

    # Check that unmapped file contains unmapped reads
    unmap_file = pysam.AlignmentFile(unmap_bam, "rb")
    unmapped_reads = [read for read in unmap_file if read.is_unmapped]
    unmap_file.close()
    assert len(unmapped_reads) > 0, "Unmapped file should contain unmapped reads"


def test_map_pe_with_unmapped_reads(tmp_path: Path):
    """Test mapping with reads that should be unmapped."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref = data_dir / "ref.fa"
    r1 = data_dir / "unmapped_read1.fq"
    r2 = data_dir / "unmapped_read2.fq"
    if not (ref.exists() and r1.exists() and r2.exists()):
        pytest.skip("unmapped test data not present")

    out_bam = str(tmp_path / "out.bam")
    unmap_bam = str(tmp_path / "unmap.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        r1_file=str(r1),
        r2_file=str(r2),
        ref_files=[str(ref)],
        output_files=out_bam,
        unmap_file=unmap_bam,
        forward_library=True,
        max_mismatches=10,
        threads=2,
        min_alignment_length=8,
        min_mapping_ratio=0.5,
        index_dir=index_dir,
        index_only=False,
        batch_size=50,
    )

    # Both files should exist
    assert Path(out_bam).exists()
    assert Path(unmap_bam).exists()

    # Check that unmapped file contains reads
    unmap_file = pysam.AlignmentFile(unmap_bam, "rb")
    unmapped_reads = list(unmap_file)
    unmap_file.close()
    assert len(unmapped_reads) > 0, "Unmapped file should contain reads"
    # Verify all reads in unmapped file are actually unmapped
    assert all(read.is_unmapped for read in unmapped_reads), (
        "All reads in unmapped file should be unmapped"
    )


def test_map_pe_multiple_refs_with_unmap_file(tmp_path: Path):
    """Test mapping with multiple references and separate unmapped file."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref1 = data_dir / "ref1.fa"
    ref2 = data_dir / "ref2.fa"
    r1 = data_dir / "test1.fq"
    r2 = data_dir / "test2.fq"
    if not all(p.exists() for p in [ref1, ref2, r1, r2]):
        pytest.skip("integration test data not present")

    out1_bam = str(tmp_path / "out1.bam")
    out2_bam = str(tmp_path / "out2.bam")
    unmap_bam = str(tmp_path / "unmap.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        r1_file=str(r1),
        r2_file=str(r2),
        ref_files=[str(ref1), str(ref2)],
        output_files=[out1_bam, out2_bam],
        unmap_file=unmap_bam,
        forward_library=True,
        max_mismatches=10,
        threads=2,
        min_alignment_length=8,
        min_mapping_ratio=0.5,
        index_dir=index_dir,
        index_only=False,
        batch_size=50,
    )

    # All files should exist
    assert Path(out1_bam).exists()
    assert Path(out2_bam).exists()
    assert Path(unmap_bam).exists()

    # Check unmapped file - verify all reads in it are unmapped (if any)
    unmap_file = pysam.AlignmentFile(unmap_bam, "rb")
    total_unmap_reads = sum(1 for _ in unmap_file)
    unmap_file.close()
    if total_unmap_reads > 0:
        # If unmapped file has reads, verify they are all unmapped
        unmap_file = pysam.AlignmentFile(unmap_bam, "rb")
        unmapped_count = sum(1 for read in unmap_file if read.is_unmapped)
        unmap_file.close()
        assert unmapped_count == total_unmap_reads, (
            "Unmapped file should only contain unmapped reads"
        )


def test_map_se_with_unmap_file(tmp_path: Path):
    """Test single-end mapping with unmapped reads written to separate file."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    ref = data_dir / "ref.fa"
    r1 = data_dir / "test1.fq"
    if not (ref.exists() and r1.exists()):
        pytest.skip("integration test data not present")

    out_bam = str(tmp_path / "out_se.bam")
    unmap_bam = str(tmp_path / "unmap_se.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        r1_file=str(r1),
        r2_file=None,
        ref_files=[str(ref)],
        output_files=out_bam,
        unmap_file=unmap_bam,
        forward_library=True,
        max_mismatches=0,
        threads=2,
        min_alignment_length=8,
        min_mapping_ratio=0.5,
        index_dir=index_dir,
        index_only=False,
        batch_size=50,
    )

    # Both files should exist
    assert Path(out_bam).exists()
    assert Path(unmap_bam).exists()

    # Check that unmapped file contains unmapped reads
    unmap_file = pysam.AlignmentFile(unmap_bam, "rb")
    unmapped_reads = [read for read in unmap_file if read.is_unmapped]
    unmap_file.close()
    assert len(unmapped_reads) > 0, "Unmapped file should contain unmapped reads"
