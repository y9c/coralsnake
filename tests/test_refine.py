"""Tests for coralsnake.refine (`coralsnake refine`)."""
import os
from pathlib import Path




FASTA = ">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 200 + "\n>chrM\n" + "G" * 50 + "\n"


def _gtf(**over):
    """A realistic 2-exon protein-coding gene with codons + UTRs."""
    return (
        'chr1\tens\tgene\t100\t300\t.\t+\t.\tgene_id "g1"; gene_name "Foo"; gene_type "protein_coding";\n'
        'chr1\tens\ttranscript\t100\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "Foo"; gene_type "protein_coding"; transcript_type "protein_coding";\n'
        'chr1\tens\texon\t100\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";\n'
        'chr1\tens\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "2";\n'
        'chr1\tens\tCDS\t110\t150\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tens\tCDS\t200\t270\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tens\tstart_codon\t110\t112\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tens\tstop_codon\t268\t270\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tens\tfive_prime_utr\t100\t109\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tens\tthree_prime_utr\t271\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )


class TestFastaRefiner:
    def test_rename_and_filter(self, tmp_path):
        from coralsnake.refine import refine_genome_references

        fa = tmp_path / "in.fa"
        fa.write_text(FASTA)
        mapper = tmp_path / "map.tsv"
        mapper.write_text("chr1\t1\nchr2\t2\nchrM\tMT\n")
        prefix = refine_genome_references(
            input_fasta=str(fa), outdir=str(tmp_path), name="test",
            rename_mapper=str(mapper),
        )
        out = Path(prefix + ".genome.fasta").read_text()
        assert out.startswith(">1 chr1")  # renamed, original kept in the description
        assert ">2 chr2" in out and ">MT chrM" in out
        # faidx + sizes produced without external samtools
        assert Path(prefix + ".genome.fasta.fai").exists()
        sizes = Path(prefix + ".genome.sizes").read_text()
        assert "1\t100" in sizes and "2\t200" in sizes

    def test_seqname_pattern_filter(self, tmp_path):
        from coralsnake.refine import refine_genome_references

        fa = tmp_path / "in.fa"
        fa.write_text(FASTA)
        prefix = refine_genome_references(
            input_fasta=str(fa), outdir=str(tmp_path), name="test",
            seqname_pattern=r"^chr[12]$",
        )
        out = Path(prefix + ".genome.fasta").read_text()
        assert ">chr1" in out and ">chr2" in out and "chrM" not in out

    def test_default_name_is_not_empty(self, tmp_path):
        """Default --outdir ./ must not produce hidden dotfile outputs."""
        from coralsnake.refine import refine_genome_references

        fa = tmp_path / "in.fa"
        fa.write_text(FASTA)
        prefix = refine_genome_references(input_fasta=str(fa), outdir=str(tmp_path))
        assert os.path.basename(prefix) != ""
        assert Path(prefix + ".genome.fasta").exists()


class TestGtfRefiner:
    def test_codon_and_utr_rows_preserved(self, tmp_path):
        """start/stop codon and UTR rows must survive refinement (metagene and
        annotate need them)."""
        from coralsnake.refine import refine_genome_references

        gtf = tmp_path / "in.gtf"
        gtf.write_text(_gtf())
        prefix = refine_genome_references(
            input_gtf=str(gtf), outdir=str(tmp_path), name="test"
        )
        out = Path(prefix + ".annotation.gtf").read_text()
        for feature in ("gene", "transcript", "exon", "CDS",
                        "start_codon", "stop_codon",
                        "five_prime_utr", "three_prime_utr"):
            assert f"\t{feature}\t" in out, f"{feature} row dropped"
        assert 'is_canonical "True"' in out

    def test_exon_only_gene_does_not_crash(self, tmp_path):
        """A gene with only exon rows gets a synthesized gene/transcript row."""
        from coralsnake.refine import refine_genome_references

        gtf = tmp_path / "in.gtf"
        gtf.write_text(
            'chr1\tens\texon\t100\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tens\texon\t200\t250\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        prefix = refine_genome_references(
            input_gtf=str(gtf), outdir=str(tmp_path), name="test"
        )
        out = Path(prefix + ".annotation.gtf").read_text()
        assert '\tgene\t' in out and '\ttranscript\t' in out

    def test_cds_without_matching_exon_is_skipped_not_crash(self, tmp_path):
        """CDS outside every exon -> transcript skipped, no crash."""
        from coralsnake.refine import refine_genome_references

        gtf = tmp_path / "in.gtf"
        gtf.write_text(
            'chr1\tens\texon\t100\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tens\tCDS\t200\t250\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        )
        prefix = refine_genome_references(
            input_gtf=str(gtf), outdir=str(tmp_path), name="test"
        )
        out = Path(prefix + ".annotation.gtf").read_text()
        skip = Path(prefix + ".skip.gtf").read_text()
        assert "g1" not in out  # the bad gene is dropped from the main GTF
        assert "g1" in skip  # ... and moved to the skip file, no crash

    def test_attribute_without_trailing_semicolon(self, tmp_path):
        """A GTF line whose last attribute lacks ';' must not crash."""
        from coralsnake.refine import refine_genome_references

        gtf = tmp_path / "in.gtf"
        gtf.write_text(
            'chr1\tens\texon\t100\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1"\n'
        )
        prefix = refine_genome_references(
            input_gtf=str(gtf), outdir=str(tmp_path), name="test"
        )
        assert "g1" in Path(prefix + ".annotation.gtf").read_text()

    def test_duplicate_gene_name_is_renamed(self, tmp_path):
        from coralsnake.refine import refine_genome_references

        gtf = tmp_path / "in.gtf"
        two = _gtf() + _gtf().replace('gene_id "g1"', 'gene_id "g2"').replace(
            'transcript_id "t1"', 'transcript_id "t2"'
        )
        gtf.write_text(two)
        prefix = refine_genome_references(
            input_gtf=str(gtf), outdir=str(tmp_path), name="test"
        )
        out = Path(prefix + ".annotation.gtf").read_text()
        assert 'gene_name "Foo_g2"' in out  # duplicate name disambiguated


class TestRefineCli:
    def test_neither_input_is_an_error(self):
        from click.testing import CliRunner

        from coralsnake.cli import cli

        res = CliRunner().invoke(cli, ["refine", "-o", "outdir"])
        assert res.exit_code != 0
        assert "Nothing to refine" in res.output
