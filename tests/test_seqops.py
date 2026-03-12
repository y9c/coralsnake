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
        assert seqops.reverse_complement("ANG") == "CNT"

    def test_empty(self):
        assert seqops.reverse_complement("") == ""

    def test_long(self):
        seq = "ACGT" * 250
        rc = seqops.reverse_complement(seq)
        assert len(rc) == 1000
        assert seqops.reverse_complement(rc) == seq


# ---------------------------------------------------------------------------
# base_conversion
# ---------------------------------------------------------------------------
class TestBaseConversion:
    def test_mk(self):
        # A→G, C→T (others unchanged)
        assert seqops.base_conversion("ACGT", "AC", "GT") == "GTGT"

    def test_km(self):
        # G→A, T→C (others unchanged)
        assert seqops.base_conversion("ACGT", "GT", "AC") == "ACAC"

    def test_lowercase(self):
        assert seqops.base_conversion("acgt", "AC", "GT") == "gtgt"

    def test_no_match(self):
        assert seqops.base_conversion("NNNN", "AC", "GT") == "NNNN"

    def test_empty(self):
        assert seqops.base_conversion("", "AC", "GT") == ""


# ---------------------------------------------------------------------------
# cal_md_and_tag
# ---------------------------------------------------------------------------
class TestCalMdAndTag:
    def test_perfect_match(self):
        md, yf, zf, yc, zc, ns, nc = seqops.cal_md_and_tag(
            "10M", "ACGTACGTAC", "ACGTACGTAC", True
        )
        assert md == "10"
        assert ns == 0

    def test_single_mismatch(self):
        md, *_ = seqops.cal_md_and_tag("4M", "ACGT", "ATGT", True)
        assert md == "1T2"

    def test_deletion(self):
        md, *_ = seqops.cal_md_and_tag("3M2D3M", "AAABBB", "AAAccBBB", True)
        assert "^" in md

    def test_insertion(self):
        md, yf, zf, yc, zc, ns, nc = seqops.cal_md_and_tag(
            "3M2I3M", "AAANNBBB", "AAABBB", True
        )
        assert md == "6"
        assert nc == 2  # insertion counted

    def test_soft_clip(self):
        md, yf, zf, yc, zc, ns, nc = seqops.cal_md_and_tag("2S3M", "NNACG", "ACG", True)
        assert md == "3"
        assert nc == 2

    def test_conversion_fwd(self):
        # fwd=True: b1=A, b2=G → ref=A, seq=G is a yf (expected conversion)
        md, yf, zf, yc, zc, ns, nc = seqops.cal_md_and_tag("2M", "GA", "AA", True)
        assert yf == 1
        assert zf == 1  # second base A matches, counted in zf

    def test_conversion_rev(self):
        # fwd=False: b1=T, b2=C → ref=T, seq=C is a yf (expected conversion)
        md, yf, zf, yc, zc, ns, nc = seqops.cal_md_and_tag("2M", "CT", "TT", False)
        assert yf == 1
        assert zf == 1

    def test_ns_count(self):
        # Non-conversion mismatch: ref=A, seq=C with fwd=True
        # b1=A, b2=G, b3=C, b4=T → seq[0]=C is neither b2 nor b4 → ns++
        md, yf, zf, yc, zc, ns, nc = seqops.cal_md_and_tag("1M", "C", "A", True)
        assert ns == 1


# ---------------------------------------------------------------------------
# calculate_directional_score
# ---------------------------------------------------------------------------
class TestDirectionalScore:
    def test_perfect(self):
        score, wrong, bad = seqops.calculate_directional_score(
            "10M", "ACGTACGTAC", "ACGTACGTAC", True
        )
        assert score == 10
        assert wrong == 0
        assert bad == 0

    def test_expected_conversion_o1(self):
        # Orientation 1: C→T is expected
        score, wrong, bad = seqops.calculate_directional_score("1M", "T", "C", True)
        assert score == 1
        assert bad == 0

    def test_wrong_conversion_o1(self):
        # Orientation 1: T→C is wrong
        score, wrong, bad = seqops.calculate_directional_score("1M", "C", "T", True)
        assert wrong == 1
        assert bad == 1

    def test_expected_conversion_o2(self):
        # Orientation 2: G→A is expected
        score, wrong, bad = seqops.calculate_directional_score("1M", "A", "G", False)
        assert score == 1
        assert bad == 0

    def test_wrong_conversion_o2(self):
        # Orientation 2: A→G is wrong
        score, wrong, bad = seqops.calculate_directional_score("1M", "G", "A", False)
        assert wrong == 1
        assert bad == 1

    def test_other_mismatch(self):
        # ref=A, seq=C, orientation 1 → other_mismatch
        score, wrong, bad = seqops.calculate_directional_score("1M", "C", "A", True)
        assert bad == 1
        assert wrong == 0

    def test_indel_penalty(self):
        # 3M + 1I + 3M = 6 matches - 1 indel = 5
        score, wrong, bad = seqops.calculate_directional_score(
            "3M1I3M", "AAANBBB", "AAABBB", True
        )
        assert score == 5

    def test_deletion_penalty(self):
        score, wrong, bad = seqops.calculate_directional_score(
            "3M2D3M", "AAABBB", "AAAccBBB", True
        )
        assert score == 4  # 6 matches - 2 deletion

    def test_soft_clip_no_penalty(self):
        score, wrong, bad = seqops.calculate_directional_score(
            "2S8M", "NNACGTACGT", "ACGTACGT", True
        )
        assert score == 8
