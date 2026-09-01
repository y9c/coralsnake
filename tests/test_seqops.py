"""Tests for the seqops C extension module."""

from coralsnake import seqops


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------
class TestReverseComplement:
    def test_palindrome(self):
        assert seqops.reverse_complement("ACGT") == "ACGT"

    def test_all_A(self):
        assert seqops.reverse_complement("AAAA") == "TTTT"

    def test_general(self):
        assert seqops.reverse_complement("ATCG") == "CGAT"

    def test_lowercase(self):
        assert seqops.reverse_complement("acgt") == "acgt"

    def test_N(self):
        # A -> T, N -> N, G -> C
        assert seqops.reverse_complement("ANG") == "CNT"

    def test_empty(self):
        assert seqops.reverse_complement("") == ""

    def test_long(self):
        seq = "ACGT" * 250
        rc = seqops.reverse_complement(seq)
        assert len(rc) == 1000
        assert seqops.reverse_complement(rc) == seq


# ---------------------------------------------------------------------------
# batch_base_conversion
# ---------------------------------------------------------------------------
class TestBatchBaseConversion:
    def test_mk(self):
        # A→G, C→T (others unchanged)
        assert seqops.batch_base_conversion(["ACGT"], "AC", "GT") == ["GTGT"]

    def test_km(self):
        # G→A, T→C (others unchanged)
        assert seqops.batch_base_conversion(["ACGT"], "GT", "AC") == ["ACAC"]

    def test_multiple(self):
        assert seqops.batch_base_conversion(["AAAA", "CCCC"], "AC", "GT") == ["GGGG", "TTTT"]


# ---------------------------------------------------------------------------
# score_and_tag
# ---------------------------------------------------------------------------
class TestScoreAndTag:
    def test_perfect_match(self):
        # score_and_tag(cigar, seq, ref, is_o1)
        # returns (score, mm, md, yf, zf, yc, zc, ns, nc)
        res = seqops.score_and_tag("10M", "ACGTACGTAC", "ACGTACGTAC", True)
        score, mm, md, yf, zf, yc, zc, ns, nc = res
        assert score == 10
        assert mm == 0
        assert md == "10"
        assert ns == 0

    def test_expected_conversion_o1(self):
        # Orientation 1: ref C -> seq T is yc (expected conversion)
        # matches: 0, exp_conv: 1, wr_conv: 0, other_mm: 0, indels: 0
        # score = matches + exp_conv - wr_conv - other_mm - indels = 0 + 1 - 0 - 0 - 0 = 1
        res = seqops.score_and_tag("1M", "T", "C", True)
        score, mm, md, yf, zf, yc, zc, ns, nc = res
        assert score == 1
        assert yc == 1
        assert mm == 0

    def test_wrong_conversion_o1(self):
        # Orientation 1: ref T -> seq C is wr_conv + ns
        # score = 0 + 0 - 1 - 0 - 0 = -1
        res = seqops.score_and_tag("1M", "C", "T", True)
        score, mm, md, yf, zf, yc, zc, ns, nc = res
        assert score == -1
        assert ns == 1
        assert mm == 1

    def test_indel_penalty(self):
        # 3M + 1I + 3M = 6 matches - 1 indel = 5
        # seq: AAANNBBB (len 8), ref: AAABBB (len 6)
        # BWA CIGAR 3M2I3M means 3 matches, 2 insertion, 3 matches
        # seq length 8 matches CIGAR 3+2+3=8
        res = seqops.score_and_tag("3M2I3M", "AAANNBBB", "AAABBB", True)
        score, mm, md, yf, zf, yc, zc, ns, nc = res
        assert score == 4 # 6 matches - 2 indels = 4

    def test_deletion_penalty(self):
        # ref: AAAccBBB (len 8), seq: AAABBB (len 6)
        # CIGAR 3M2D3M
        res = seqops.score_and_tag("3M2D3M", "AAABBB", "AAAccBBB", True)
        score, mm, md, yf, zf, yc, zc, ns, nc = res
        assert score == 4 # 6 matches - 2 indels = 4
        assert md == "3^cc3"


# reverse_md (C kernel added for tbam2gbam liftover)
class TestReverseMd:
    # reference vectors mirror tbam2gbam.test_tbam2gbam / the Python implementation
    def test_known_cases(self):
        cases = {
            "10": "10",
            "0": "0",
            "2A3": "3T2",
            "0A0C5": "5G0T0",
            "2^ACG3": "3^CGT2",
            "0A0^GT3": "3^AC0T0",
            "0G0": "0C0",
            "3A2^TT0C4": "4G0^AA2T3",
            "0^ACGT0": "0^ACGT0",
        }
        for md, expected in cases.items():
            assert seqops.reverse_md(md) == expected, md

    def test_matches_python_reference(self):
        import re

        COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

        def py_reverse_md(md):
            parts = []
            for m in re.finditer(r"([0-9]+)|([A-Z])|(\^[A-Z]+)", md):
                if m.group(1) is not None:
                    parts.append(int(m.group(1)))
                elif m.group(2):
                    parts.append(m.group(2).translate(COMP))
                elif m.group(3):
                    del_bases = m.group(3)[1:]
                    parts.append("^" + del_bases.translate(COMP)[::-1])
            parts.reverse()
            merged, cur = [], 0
            for p in parts:
                if isinstance(p, int):
                    cur += p
                else:
                    merged.append(str(cur))
                    cur = 0
                    merged.append(p)
            merged.append(str(cur))
            return "".join(merged)

        for md in [
            "0", "5", "100", "1A0", "47T13", "5N3", "10^ACGT2", "3G2C1",
            "0A^GT0", "2^C3TA", "1^AC2G3", "7A^CT6^G5", "0^NN0",
        ]:
            assert seqops.reverse_md(md) == py_reverse_md(md), md
