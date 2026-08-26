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
        tx = Transcript(
            exons={1: __import__("coralsnake.utils", fromlist=["Span"]).Span(0, 500)},
            strand="+",
        )
        tx.priority = (10, 0)
        r = rank_transcript("UNKNOWN", tx)
        assert r[0] == 4


# ---------------------------------------------------------------------------
# Transcript.to_tsv with_txpos
# ---------------------------------------------------------------------------
class TestToTsvTxPos:
    def test_plus_bounding_box(self):
        from coralsnake.utils import Transcript, Span

        tx = Transcript(gene_id="G", transcript_id="T", chrom="c")
        tx.add_exon("1", Span(10, 20))
        tx.add_exon("2", Span(30, 55))
        s = tx.to_tsv(with_txpos=True)
        fields = s.split("\t")
        # spans, start_codon/stop_codon skipped (with_codon=False) => last 2 = txpos
        assert fields[-2:] == ["11", "55"]  # min start+1, max end

    def test_minus_bounding_box_is_valid_span(self):
        """Regression: gene on '-' strand used to emit start>end (inverted).

        transcript_start/transcript_end are the genomic bounding box
        [min_start, max_end], identical for both strands.
        """
        from coralsnake.utils import Transcript, Span

        # two exons, 5'->3' = right to left; unsorted insertion order
        tx = Transcript(gene_id="G", transcript_id="T", chrom="c", strand="-")
        tx.add_exon("1", Span(151096, 151166))  # 5' exon (rightmost)
        tx.add_exon("2", Span(147593, 151006))  # 3' exon (leftmost)
        s = tx.to_tsv(with_txpos=True)
        fields = s.split("\t")
        start, end = int(fields[-2]), int(fields[-1])
        assert start <= end
        assert (start, end) == (147594, 151166)  # min start+1, max end


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
