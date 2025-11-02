from pathlib import Path

import pytest
import pysam

from coralsnake.tbam2gbam import convert_bam


def test_liftover_basic(tmp_path: Path):
    """Test liftover command - convert transcript BAM to genomic BAM."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    input_bam = data_dir / "liftover_input.bam"
    annot_file = data_dir / "liftover_annot.tsv"
    faidx_file = data_dir / "liftover_faidx.fai"

    if not all(p.exists() for p in [input_bam, annot_file, faidx_file]):
        pytest.skip("liftover test data not present")

    output_bam = str(tmp_path / "output.bam")

    convert_bam(
        input_bam=str(input_bam),
        output_bam=output_bam,
        annotation_file=str(annot_file),
        faidx_file=str(faidx_file),
        threads=2,
        sort=False,
    )

    # Check output file exists
    assert Path(output_bam).exists()

    # Check output BAM is readable
    with pysam.AlignmentFile(output_bam, "rb") as f:
        # Should be able to read the header at least
        assert f.header is not None
