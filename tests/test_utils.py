"""Tests for coralsnake.utils – Transcript, Span, helpers."""

from coralsnake.utils import (
    Span,
    Transcript,
    format_duration,
    load_annotation,
    load_faidx,
)


# ---------------------------------------------------------------------------
# Span
# ---------------------------------------------------------------------------
class TestSpan:
    def test_fields(self):
        s = Span(10, 20)
        assert s.start == 10
        assert s.end == 20

    def test_iter(self):
        assert list(Span(5, 15)) == [5, 15]

    def test_repr(self):
        assert "start=5" in repr(Span(5, 15))


# ---------------------------------------------------------------------------
# Transcript – construction and properties
# ---------------------------------------------------------------------------
class TestTranscript:
    def _make(self, strand="+"):
        return Transcript(
            gene_id="G1",
            transcript_id="T1",
            chrom="chr1",
            strand=strand,
            exons={
                1: Span(100, 200),  # 100 bp
                2: Span(300, 450),  # 150 bp
                3: Span(500, 530),  #  30 bp
            },
        )

    # --- exons_forwards ---
    def test_exons_forwards_plus(self):
        tx = self._make("+")
        fwd = tx.exons_forwards
        assert [e.start for e in fwd] == [100, 300, 500]

    def test_exons_forwards_minus(self):
        tx = self._make("-")
        fwd = tx.exons_forwards
        assert [e.start for e in fwd] == [100, 300, 500]

    # --- cum_exon_lens ---
    def test_cum_exon_lens_plus(self):
        assert self._make("+").cum_exon_lens == [100, 250, 280]

    def test_cum_exon_lens_minus(self):
        assert self._make("-").cum_exon_lens == [100, 250, 280]

    def test_cum_exon_lens_cached(self):
        tx = self._make()
        a = tx.cum_exon_lens
        b = tx.cum_exon_lens
        assert a is b

    def test_cum_exon_lens_invalidated_on_add(self):
        tx = self._make()
        _ = tx.cum_exon_lens
        tx.add_exon(4, Span(600, 650))
        assert tx.cum_exon_lens[-1] == 330

    def test_cum_exon_lens_invalidated_on_sort(self):
        tx = self._make()
        _ = tx.cum_exon_lens
        tx.sort_exons()
        assert tx._cum_exon_lens is None  # cleared by sort

    # --- length ---
    def test_length(self):
        assert self._make().length == 280

    def test_length_single_exon(self):
        tx = Transcript(exons={1: Span(0, 50)}, strand="+")
        assert tx.length == 50

    # --- to_tsv / get_genome_spans ---
    def test_to_tsv_basic(self):
        tx = self._make()
        tsv = tx.to_tsv()
        fields = tsv.split("\t")
        assert fields[0] == "G1"
        assert fields[1] == "T1"
        assert "101-200" in fields[4]

    def test_get_genome_spans(self):
        tx = self._make()
        spans = tx.get_genome_spans()
        assert "101-200" in spans

    def test_get_gene_spans(self):
        tx = self._make()
        spans = tx.get_gene_spans()
        assert spans.startswith("1-100")

    # --- repr ---
    def test_repr(self):
        tx = self._make()
        r = repr(tx)
        assert "G1" in r and "T1" in r


# ---------------------------------------------------------------------------
# load_annotation / load_faidx
# ---------------------------------------------------------------------------
class TestLoadAnnotation:
    def test_round_trip(self, tmp_path):
        tsv = tmp_path / "annot.tsv"
        tsv.write_text(
            "gene_id\ttranscript_id\tchrom\tstrand\tspans\n"
            "G1\tT1\tchr1\t+\t101-200,301-400\n"
            "G1\tT2\tchr1\t-\t501-600\n"
        )
        annot = load_annotation(str(tsv))
        assert "G1" in annot
        assert "T1" in annot["G1"]
        assert "T2" in annot["G1"]
        t1 = annot["G1"]["T1"]
        assert t1.chrom == "chr1"
        assert t1.strand == "+"
        assert len(t1.exons) == 2

    def test_no_header(self, tmp_path):
        tsv = tmp_path / "annot.tsv"
        tsv.write_text("G1\tT1\tchr1\t+\t101-200\n")
        annot = load_annotation(str(tsv), with_header=False)
        assert "T1" in annot["G1"]


class TestLoadFaidx:
    def test_basic(self, tmp_path):
        fai = tmp_path / "ref.fai"
        fai.write_text("chr1\t1000\t5\t80\t81\nchr2\t2000\t1100\t80\t81\n")
        idx = load_faidx(str(fai))
        assert idx == {"chr1": 1000, "chr2": 2000}


# ---------------------------------------------------------------------------
# format_duration
# ---------------------------------------------------------------------------
class TestFormatDuration:
    def test_seconds(self):
        assert format_duration(5.123) == "5.12s"

    def test_minutes(self):
        d = format_duration(125)
        assert "2m" in d

    def test_hours(self):
        d = format_duration(3700)
        assert "1h" in d
