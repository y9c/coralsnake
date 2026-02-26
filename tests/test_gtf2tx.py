"""Tests for coralsnake.gtf2tx – GTF/GFF parsing and transcript ranking."""

from coralsnake.gtf2tx import (
    parse_gff_annot,
    parse_gtf_annot,
    rank_transcript,
    sanitize_sequence_name,
)
from coralsnake.utils import Transcript


# ---------------------------------------------------------------------------
# GTF / GFF attribute parsing
# ---------------------------------------------------------------------------
class TestParseGtfAnnot:
    def test_basic(self):
        d = parse_gtf_annot('gene_id "ENSG001"; transcript_id "ENST001";')
        assert d["gene_id"] == "ENSG001"
        assert d["transcript_id"] == "ENST001"

    def test_duplicate_keys(self):
        d = parse_gtf_annot('tag "basic"; tag "MANE_Select";')
        assert "basic" in d["tag"]
        assert "MANE_Select" in d["tag"]

    def test_empty(self):
        d = parse_gtf_annot("")
        assert d == {}


class TestParseGffAnnot:
    def test_basic(self):
        d = parse_gff_annot("ID=exon1;Parent=mRNA1")
        assert d["ID"] == "exon1"
        assert d["Parent"] == "mRNA1"

    def test_empty(self):
        d = parse_gff_annot("")
        assert d == {}


# ---------------------------------------------------------------------------
# rank_transcript
# ---------------------------------------------------------------------------
class TestRankTranscript:
    def test_mane_select_highest(self):
        tx = Transcript()
        tx.priority = (0, 0)
        assert rank_transcript("T1", tx) == (0, 0)

    def test_ensembl_canonical(self):
        tx = Transcript()
        tx.priority = (0, 1)
        assert rank_transcript("T1", tx) == (0, 1)

    def test_dot_suffix(self):
        tx = Transcript()
        tx.priority = (10, 0)
        assert rank_transcript("AT1G01010.1", tx) == (2, 1)

    def test_dash_suffix(self):
        tx = Transcript()
        tx.priority = (10, 0)
        assert rank_transcript("LOC_Os01g01010-01", tx) == (3, 1)

    def test_fallback_length(self):
        tx = Transcript(exons={1: __import__("coralsnake.utils", fromlist=["Span"]).Span(0, 500)}, strand="+")
        tx.priority = (10, 0)
        r = rank_transcript("UNKNOWN", tx)
        assert r[0] == 4


# ---------------------------------------------------------------------------
# sanitize_sequence_name
# ---------------------------------------------------------------------------
class TestSanitizeSequenceName:
    def test_clean_name(self):
        assert sanitize_sequence_name("ENSG00000001") == "ENSG00000001"

    def test_special_chars(self):
        result = sanitize_sequence_name("gene (copy)")
        assert "(" not in result
        assert "_" in result

    def test_empty(self):
        assert sanitize_sequence_name("") == ""


# ---------------------------------------------------------------------------
# Integration: parse_file
# ---------------------------------------------------------------------------
class TestParseFile:
    def test_basic(self, tmp_path, data_dir, has_gtf_data):
        from pathlib import Path
        from coralsnake.gtf2tx import parse_file

        output = str(tmp_path / "output.tsv")
        parse_file(
            gtf_file=str(data_dir / "R64-1-1.release57.gtf"),
            fasta_file=str(data_dir / "R64-1-1.fa"),
            output_file=output,
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
        content = Path(output).read_text()
        lines = content.strip().split("\n")
        assert len(lines) > 1  # header + data

    def test_with_seq(self, tmp_path, data_dir, has_gtf_data):
        from pathlib import Path
        from coralsnake.gtf2tx import parse_file

        output = str(tmp_path / "output.tsv")
        seq_out = str(tmp_path / "seqs.fa")
        parse_file(
            gtf_file=str(data_dir / "R64-1-1.release57.gtf"),
            fasta_file=str(data_dir / "R64-1-1.fa"),
            output_file=output,
            seq_file=seq_out,
            sanitize=False,
            with_codon=False,
            with_genename=False,
            with_biotype=False,
            with_txpos=False,
            filter_biotype=None,
            seq_upper=True,
            line_length=0,
        )
        assert Path(seq_out).exists()
        seq_content = Path(seq_out).read_text()
        assert seq_content.startswith(">")
