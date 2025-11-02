from pathlib import Path

import pytest

from coralsnake.genegroup import group_genes


def test_group_basic(tmp_path: Path):
    """Test group command - group and find consensus of genes."""
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    gtf_file = data_dir / "R64-1-1.release57.gtf"
    fasta_file = data_dir / "R64-1-1.fa"

    if not (gtf_file.exists() and fasta_file.exists()):
        pytest.skip("group test data not present")

    output_file = str(tmp_path / "output.tsv") if tmp_path else None
    output_consensus = str(tmp_path / "consensus.fa") if tmp_path else None

    group_genes(
        fa_file_list=[str(fasta_file)],
        gtf_file_list=[str(gtf_file)],
        out_file=output_file,
        consensus_fa=output_consensus,
        gene_name_regex=None,
        gene_biotype_list=None,
        gene_length_limit=300,
        cluster_threshold=0.1,
        threads=2,
    )

    # If output file specified, check it exists
    if output_file:
        assert Path(output_file).exists()
        content = Path(output_file).read_text()
        assert len(content) > 0

    # If consensus file specified, check it exists
    if output_consensus:
        assert Path(output_consensus).exists()
        content = Path(output_consensus).read_text()
        assert len(content) > 0
