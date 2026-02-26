"""Tests for coralsnake.genegroup – rename helpers and reassignment logic."""

import re

import pytest

from coralsnake.genegroup import rename_snRNA, rename_snoRNA


# ---------------------------------------------------------------------------
# rename_snRNA
# ---------------------------------------------------------------------------
class TestRenameSnRNA:
    def test_RNU2(self):
        assert rename_snRNA("RNU2-1") == "U2"

    def test_RNU2_pseudo(self):
        assert rename_snRNA("RNU2-27P") == "U2"

    def test_RNU4ATAC(self):
        assert rename_snRNA("RNU4ATAC7") == "U4atac"

    def test_RNVU1(self):
        assert rename_snRNA("RNVU1-1") == "U1V"

    def test_drosophila(self):
        assert rename_snRNA("snRNA:U2:38ABb") == "U2"

    def test_passthrough(self):
        assert rename_snRNA("unrelated") == "unrelated"


# ---------------------------------------------------------------------------
# rename_snoRNA
# ---------------------------------------------------------------------------
class TestRenameSnoRNA:
    def test_ACA(self):
        assert rename_snoRNA("ACA64") == "SNORA64"

    def test_U8(self):
        assert rename_snoRNA("U8") == "SNORD118"

    def test_U3(self):
        r = rename_snoRNA("U3")
        assert r.startswith("SNORD")

    def test_Z6(self):
        assert rename_snoRNA("Z6") == "Z6"


# ---------------------------------------------------------------------------
# Reassignment logic (collect-then-apply)
# ---------------------------------------------------------------------------
class TestReassignment:
    """Verify the collect-then-apply pattern doesn't skip items."""

    def test_no_skip(self):
        class FakeTx:
            def __init__(self, name, seq):
                self.name = name
                self._seq = seq

        named = {
            "geneA": [FakeTx("a1", "ACGTACGTACGT")],
            "geneB": [FakeTx("b1", "GGGGGGGG")],
            "unnamed": [
                FakeTx("u1", "ACGT"),       # too short ratio for geneA
                FakeTx("u2", "ACGTACGT"),    # matches geneA (8 > 0.5*12)
                FakeTx("u3", "GGGG"),        # too short ratio for geneB
            ],
        }

        reassigned = []
        for tx in named["unnamed"]:
            matched = False
            for target_gene, tx_list in named.items():
                if target_gene == "unnamed":
                    continue
                for tx2 in tx_list:
                    if tx._seq in tx2._seq and len(tx._seq) > 0.5 * len(tx2._seq):
                        reassigned.append((tx, target_gene))
                        matched = True
                        break
                if matched:
                    break

        for tx, target_gene in reassigned:
            named["unnamed"].remove(tx)
            named[target_gene].append(tx)

        assert len(named["unnamed"]) == 2
        assert len(named["geneA"]) == 2
        assert any(tx.name == "u2" for tx in named["geneA"])


# ---------------------------------------------------------------------------
# Integration: group_genes (requires test data)
# ---------------------------------------------------------------------------
class TestGroupGenes:
    def test_basic(self, tmp_path, data_dir, has_gtf_data):
        from coralsnake.genegroup import group_genes

        output = str(tmp_path / "output.tsv")
        consensus = str(tmp_path / "consensus.fa")
        group_genes(
            fa_file_list=[str(data_dir / "R64-1-1.fa")],
            gtf_file_list=[str(data_dir / "R64-1-1.release57.gtf")],
            out_file=output,
            consensus_fa=consensus,
            threads=2,
        )
        from pathlib import Path
        assert Path(output).exists()
        assert Path(output).stat().st_size > 0
