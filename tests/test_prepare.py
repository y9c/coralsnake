from pathlib import Path

import pytest

from coralsnake.gtf2tx import parse_file


def test_prepare_basic(tmp_path: Path):
    """Test prepare command - extract primary transcript from GTF."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    gtf_file = data_dir / "R64-1-1.release57.gtf"
    fasta_file = data_dir / "R64-1-1.fa"

    if not (gtf_file.exists() and fasta_file.exists()):
        pytest.skip("prepare test data not present")

    output_file = str(tmp_path / "output.tsv")

    parse_file(
        gtf_file=str(gtf_file),
        fasta_file=str(fasta_file),
        output_file=output_file,
        seq_file=None,
        sanitize=False,
        with_codon=False,
        with_genename=False,
        with_biotype=False,
        with_txpos=False,
        filter_biotype=None,
        seq_upper=True,
        line_length=0,
    )

    # Check output file exists and has content
    assert Path(output_file).exists()
    content = Path(output_file).read_text()
    assert len(content) > 0
    # Should have header and at least one transcript
    lines = content.strip().split("\n")
    assert len(lines) > 1
