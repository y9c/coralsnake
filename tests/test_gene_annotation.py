"""Tests for coralsnake.gene_annotation — the shared gene-annotation model."""

import os

from coralsnake.gene_annotation import (
    AnnotationRow,
    GeneAnnotation,
    parse_gff_annot,
    parse_gtf_annot,
)


GTF = (
    "##gtf-version 2.2\n"
    "#!genome-build R64-1-1\n"
    'chr1\tens\tgene\t100\t500\t.\t+\t.\tgene_id "g1"; gene_name "Foo"; gene_biotype "protein_coding";\n'
    'chr1\tens\ttranscript\t100\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
    'chr1\tens\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";\n'
    'chr1\tens\texon\t300\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "2";\n'
    'chr1\tens\tCDS\t110\t200\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tens\tstart_codon\t110\t112\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tens\tstop_codon\t488\t490\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    'chr2\tens\tgene\t10\t90\t.\t-\t.\tgene_id "g2"; gene_name "Bar";\n'
    'chr2\tens\ttranscript\t10\t90\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'
    'chr2\tens\texon\t10\t90\t.\t-\t.\tgene_id "g2"; transcript_id "t2"; exon_number "1";\n'
)


def _write(tmp_path, text, name="in.gtf"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


class TestAttributeParsers:
    def test_parsers_reexported_from_gtf2tx(self):
        from coralsnake.gtf2tx import parse_gff_annot as g2
        from coralsnake.gtf2tx import parse_gtf_annot as g1

        assert g1('gene_id "G1"; transcript_id "T1";') == parse_gtf_annot(
            'gene_id "G1"; transcript_id "T1";'
        )
        assert g2("ID=exon1;Parent=mRNA1") == parse_gff_annot("ID=exon1;Parent=mRNA1")

    def test_missing_trailing_semicolon_and_quoted_semicolon(self):
        d = parse_gtf_annot('gene_id "G1"; gene "a; b"')
        assert d == {"gene_id": "G1", "gene": "a; b"}


class TestGeneAnnotationParse:
    def test_grouping_and_accessors(self, tmp_path):
        ann = GeneAnnotation(_write(tmp_path, GTF))
        assert list(ann.genes) == ["g1", "g2"]
        assert list(ann.iter_rows())[0].feature == "gene"

        g1 = ann.gene("g1")
        assert g1.seqname == "chr1" and g1.strand == "+"
        assert g1.span == (100, 500)
        assert len(g1.transcripts) == 1

        tx = g1.transcripts["t1"]
        assert len(tx.exons) == 2 and len(tx.cds) == 1
        assert tx.is_coding and not ann.gene("g2").transcripts["t2"].is_coding
        assert len(tx.start_codons) == 1 and len(tx.stop_codons) == 1
        # 0-based half-open spans
        assert tx.exons[0].span_0 == (99, 200)
        assert tx.exons[0].length == 101

    def test_attributes_are_lossless(self, tmp_path):
        ann = GeneAnnotation(_write(tmp_path, GTF))
        row = ann.gene("g1").gene_rows[0]
        assert row.attributes["gene_name"] == "Foo"

    def test_rows_without_gene_id_are_kept_unassigned(self, tmp_path):
        text = GTF + "chr1\tfoo\texon\t5\t9\t.\t+\t.\tno_gene here\n"
        ann = GeneAnnotation(_write(tmp_path, text))
        assert len(ann.unassigned_rows) == 1
        # 10 GTF data rows + 1 unassigned row (header lines are not rows)
        assert len(list(ann.iter_rows())) == 11
        # pruning applies the GTF-cleaning policy
        removed = ann.prune()
        assert removed == 1
        assert len(ann.unassigned_rows) == 0

    def test_malformed_lines_skipped(self, tmp_path):
        text = "chr1\tfoo\n" + GTF
        ann = GeneAnnotation(_write(tmp_path, text))
        assert list(ann.genes) == ["g1", "g2"]

    def test_gzip_input(self, tmp_path):
        import gzip

        p = tmp_path / "in.gtf.gz"
        with gzip.open(p, "wt") as f:
            f.write(GTF)
        ann = GeneAnnotation(str(p))
        assert list(ann.genes) == ["g1", "g2"]

    def test_seqname_mapper_and_pattern(self, tmp_path):
        ann = GeneAnnotation(
            _write(tmp_path, GTF),
            seqname_mapper={"chr1": "1", "chr2": "2"},
            seqname_pattern=r"^[12]$",
        )
        assert list(ann.genes) == ["g1", "g2"]
        assert ann.gene("g1").seqname == "1"
        ann2 = GeneAnnotation(_write(tmp_path, GTF), seqname_pattern=r"^chr1$")
        assert list(ann2.genes) == ["g1"]

    def test_gff_attribute_style(self, tmp_path):
        gff = (
            "##gff-version 3\n"
            "chr1\tsgd\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=Foo\n"
            "chr1\tsgd\tmRNA\t100\t200\t.\t+\t.\tID=mrna1;Parent=gene1\n"
        )
        ann = GeneAnnotation(_write(tmp_path, gff, "in.gff3"))
        assert ann.is_gff
        # GFF3 rows carry ID/Parent (assigned via parse_gff_annot), no gene_id
        assert "gene1" not in ann.genes
        assert len(ann.unassigned_rows) == 2


class TestGeneAnnotationWrite:
    def test_roundtrip(self, tmp_path):
        ann = GeneAnnotation(_write(tmp_path, GTF))
        out = str(tmp_path / "out.gtf")
        ann.write_gtf(out, sort=False, bgzip=False)
        txt = open(out).read()
        # header preserved
        assert txt.startswith("##gtf-version 2.2\n")
        # attributes preserved verbatim
        assert 'gene_name "Foo"' in txt and 'exon_number "1"' in txt
        # body in original file order (gene1 before gene2)
        assert txt.index('gene_id "g1"') < txt.index('gene_id "g2"')
        # re-parsing yields the same model
        ann2 = GeneAnnotation(out)
        assert list(ann2.genes) == ["g1", "g2"]
        assert ann2.gene("g1").transcripts["t1"].cds[0].start == 110

    def test_sort_and_bgzip(self, tmp_path):
        # write g2 before g1, then sort puts chr1 gene first
        ann = GeneAnnotation(_write(tmp_path, GTF))
        g2 = ann.genes.pop("g2")
        ann.genes = {"g2": g2, **ann.genes}
        out = str(tmp_path / "out.gtf")
        ann.write_gtf(out, sort=True, bgzip=True)
        lines = [ln for ln in open(out) if not ln.startswith("#")]
        assert lines[0].startswith("chr1\t")
        assert len(lines) == 10  # 7 g1 rows + 3 g2 rows
        assert any(ln.startswith("chr2\t") for ln in lines)
        assert open(out + ".gz", "rb").read(2) == b"\x1f\x8b"  # gzip magic
        assert os.path.exists(out + ".gz.tbi")

    def test_remove_gene(self, tmp_path):
        ann = GeneAnnotation(_write(tmp_path, GTF))
        ann.remove_gene("g2")
        assert list(ann.genes) == ["g1"]
        assert all(r.attributes.get("gene_id") != "g2" for r in ann.iter_rows())


class TestAnnotationRow:
    def test_serialize(self):
        row = AnnotationRow(
            "chr1", "ens", "exon", 100, 200, ".", "+", ".", {"gene_id": "g"}
        )
        assert row.to_gtf_line() == 'chr1\tens\texon\t100\t200\t.\t+\t.\tgene_id "g";'
        assert row.span_0 == (99, 200)
        assert row.length == 101
