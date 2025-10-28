from pathlib import Path

import pytest
import pysam

from coralsnake.mapping import map_file


def test_map_pe_integration(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests"
    ref = data_dir / "ref.fa"
    r1 = data_dir / "test1.fq"
    r2 = data_dir / "test2.fq"
    if not (ref.exists() and r1.exists() and r2.exists()):
        pytest.skip("integration test data not present")

    out_bam = str(tmp_path / "out.bam")
    index_dir = str(tmp_path / "idx")

    map_file(
        ref_file=str(ref),
        r1_file=str(r1),
        r2_file=str(r2),
        output_file=out_bam,
        fwd_lib=True,
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
    data_dir = repo_root / "tests"
    ref = data_dir / "ref.fa"
    if not ref.exists():
        pytest.skip("integration test data not present")

    index_dir = tmp_path / "idx"
    map_file(
        ref_file=str(ref),
        r1_file=None,
        r2_file=None,
        output_file=str(tmp_path / "dummy.bam"),
        index_dir=str(index_dir),
        index_only=True,
    )

    # BWA index should exist (ref.orig.fa + ref.mk.fa*.{amb,ann,bwt,pac,sa})
    assert (index_dir / "ref.orig.fa").exists()
    assert (index_dir / "ref.mk.amb").exists()
