from pathlib import Path

from coralsnake.annot import run_annot


def write_file(p: Path, content: str) -> str:
    p.write_text(content)
    return str(p)


def test_run_annot_basic(tmp_path: Path):
    # Build a tiny annotation: two exons on chr1 + strand: 10-20 and 30-40 (1-based in file)
    annot_path = write_file(
        tmp_path / "annot.tsv",
        "\t".join(["chrom", "strand", "spans", "gene_id", "transcript_id"])
        + "\n"
        + "\t".join(["chr1", "+", "10-20,30-40", "GENE1", "TX1"])
        + "\n",
    )

    # Input sites: columns chrom pos strand (pos 1-based)
    # Expect: chr1 10 + => transcript_pos 0 (first base of first exon)
    #         chr1 35 + => transcript_pos 16 (within second exon; first exon length 11)
    #         chr1 5  + => NA
    input_path = write_file(
        tmp_path / "sites.tsv",
        "\n".join(
            [
                "chr1\t10\t+",
                "chr1\t35\t+",
                "chr1\t5\t+",
            ]
        )
        + "\n",
    )

    out_path = str(tmp_path / "out.tsv")

    run_annot(
        input_file=input_path,
        output_file=out_path,
        annot_file=annot_path,
        cols=None,
        keep_na=True,
        collapse_annot=False,
        add_count=False,
        skip_header=False,
    )

    lines = Path(out_path).read_text().rstrip().splitlines()
    assert lines[0] == "chr1\t10\t+\tGENE1\tTX1\t0"
    assert lines[1] == "chr1\t35\t+\tGENE1\tTX1\t16"
    assert lines[2] == "chr1\t5\t+\t.\t.\t."
