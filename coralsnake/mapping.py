#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-08 20:32


import os
import random

from . import mappy as mp
from .conversion import km_conversion, mk_conversion, mk_convert_file


def find_properly_paired_hits(hits, fwd=True):
    """
    # find properly paired hits
    Given a list of hits with read_num, ref_name, ref_start, ref_end, strand. For example:
        1 rRNA-Hsa-nucleus_locus 8204 8237 1
        1 rRNA-Hsa-nucleus_locus 21555 21588 1
        2 rRNA-Hsa-nucleus_locus 21561 21596 -1
        2 rRNA-Hsa-nucleus_locus 8210 8245 -1

    This function will find all properly paired hits. A hit is properly paired if it is properly paired.
    We defined properly paired hits as:
    1. The hit paris are from two reads (read_num=1 and read_num=2)
    2. Both hits mapped to the same reference sequence
    3. read 1 and reads 2 of one pair are in different strand (if reads 1 is 1, reads 2 is -1)
    4. If fwd is True, strand = 1 is forward, strand = -1 is reverse. Otherwise, strand = 1 is reverse, strand = -1 is forward
    5. The mapping start of read 1 is smaller than the mapping end of read 2, if fwd is True. Otherwise, the mapping end of read 1 is larger than the mapping start of read 2.
    6. The distance between the mapping start and end of read 1 and read 2 is less than 1000
    """
    parsed_hits = []
    # group by ref_name and separate read 1 and read 2
    ref_name_hits = {}
    for hit in hits:
        if hit.ctg not in ref_name_hits:
            ref_name_hits[hit.ctg] = [[], []]
        ref_name_hits[hit.ctg][hit.read_num - 1].append(hit)
    for hits in ref_name_hits.values():
        if len(hits[0]) > 0 and len(hits[1]) > 0:
            for hit1 in hits[0]:
                for hit2 in hits[1]:
                    if hit1.strand + hit2.strand == 0:
                        if fwd:
                            if hit1.r_st < hit2.r_en and hit2.r_en - hit1.r_st < 1000:
                                parsed_hits.append((hit1, hit2))
                        else:
                            if hit1.r_en > hit2.r_st and hit1.r_en - hit2.r_st < 1000:
                                parsed_hits.append((hit1, hit2))

    return parsed_hits


def cal_md_and_tag(cigar, seq, ref, fwd):
    """
    fwd:
    if True, A->G, C->T
    if False, T->C, G->A

    calculate MD tag
    M (0): Alignment match
    I (1): Insertion to the reference
    D (2): Deletion from the reference
    N (3): Skipped region from the reference
    S (4): Soft clipping
    H (5): Hard clipping
    P (6): Padding
    = (7): Sequence match
    X (8): Sequence mismatch

    Some tags:
    Yf:i:<N>: Number of (A to G) conversions are detected in the read.
    Zf:i:<N>: Number of (A to A) un-converted bases are detected in the read. Yf + Zf = total number of bases which can be converted in the read sequence.
    Yc:i:<N>: Number of (C to T) conversions are detected in the read.
    Zc:i:<N>: Number of (C to C) un-converted bases are detected in the read. Yc + Zc = total number of bases which can be converted in the read
    NS:i:<N>: Number of substitutions are detected in the read.
    NC:i:<N>: Number of clipped bases and INDEL bases are detected in the read.


    ### add tag function that I used previously, learn from this function
    for read_pos, _, ref_base in read.get_aligned_pairs(with_seq=True):
        # Dn not forget to convert the ref_base to upper case
        ref_base = ref_base.upper() if ref_base is not None else None
        read_base = s[read_pos] if read_pos is not None else None
        if strand == "+":
            if ref_base == "A":
                if read_base == "G":
                    yf += 1
                elif read_base == "A":
                    zf += 1
                elif read_base is not None:
                    ns += 1
            elif ref_base == "C":
                if read_base == "T":
                    yc += 1
                elif read_base == "C":
                    zc += 1
                elif read_base is not None:
                    ns += 1
            elif ref_base is None or read_base is None:
                nc += 1
            elif ref_base != read_base:
                ns += 1

        else:
            if ref_base == "T":
                if read_base == "C":
                    yf += 1
                elif read_base == "T":
                    zf += 1
                elif read_base is not None:
                    ns += 1
            elif ref_base == "G":
                if read_base == "A":
                    yc += 1
                elif read_base == "G":
                    zc += 1
                elif read_base is not None:
                    ns += 1
            elif ref_base is None or read_base is None:
                nc += 1
            elif ref_base != read_base:
                ns += 1
    """
    yf = 0
    zf = 0
    yc = 0
    zc = 0
    ns = 0
    nc = 0
    md_tag = []
    ref_index = 0
    query_index = 0
    match_count = 0

    if fwd:
        b1, b2, b3, b4 = "A", "G", "C", "T"
    else:
        b1, b2, b3, b4 = "T", "C", "G", "A"

    for length, operation in cigar:
        if operation == 0:  # Match or Mismatch
            for i in range(length):
                if ref[ref_index] == seq[query_index]:
                    match_count += 1
                    if seq[query_index] == b1:
                        zf += 1
                    elif seq[query_index] == b3:
                        zc += 1
                else:
                    if match_count > 0:
                        md_tag.append(str(match_count))
                        match_count = 0
                    md_tag.append(ref[ref_index])
                    if seq[query_index] == b2:
                        yf += 1
                    elif seq[query_index] == b4:
                        yc += 1
                    else:
                        ns += 1
                ref_index += 1
                query_index += 1
        elif operation == 1:  # Insertion to the reference (ignored in MD tag)
            query_index += length
            nc += length
        elif operation == 4:  # Soft clipping
            query_index += length
            nc += length
        elif operation == 2:  # Deletion from the reference
            if match_count > 0:
                md_tag.append(str(match_count))
                match_count = 0
            md_tag.append("^" + ref[ref_index : ref_index + length])
            ref_index += length
            nc += length

    if match_count > 0:
        md_tag.append(str(match_count))

    return "".join(md_tag), yf, zf, yc, zc, ns, nc


def filter_hits(hits, seq1, seq2):
    """
    Filter hits by the following rules:
    1. The mapping quality of the hit is greater than 0
    2. The alignment length of the hit is greater than 20
    3. The mapping length of the hit is larger than 80% of the query length
    """
    filtered_hits = []
    for hit in hits:
        q_len = len(seq1) if hit.read_num == 1 else len(seq2)
        if hit.mapq > 0 and hit.blen > 20 and hit.mlen > 0.8 * q_len:
            filtered_hits.append(hit)
    return filtered_hits


def run_mapping(name, seq1, seq2, qua1, qua2, idx0, idx_km, idx_mk, fwd_lib=True):
    # sam format specification
    # 1. QNAME: Query template NAME
    # 2. FLAG: bitwise FLAG
    # 3. RNAME: Reference sequence NAME
    # 4. POS: 1-based leftmost POSition/coordinate of clipped sequence
    # 5. MAPQ: MAPping Quality (Phred-scaled)
    # 6. CIGAR: CIGAR string
    # 7. RNEXT: Reference name of the mate/next read
    # 8. PNEXT: Position of the mate/next read
    # 9. TLEN: observed Template LENgth
    # 10. SEQ: segment SEQuence
    # 11. QUAL: ASCII of Phred-scaled base QUALity+33
    # 12. TAG: additional information

    # q_st  q_en  strand  ctg  ctg_len  r_st  r_en  mlen  blen  mapq  cg:Z:cigar_str
    if fwd_lib:
        seq1_conv = mk_conversion(seq1)
        seq2_conv = km_conversion(seq2)
    else:
        seq1_conv = km_conversion(seq1)
        seq2_conv = mk_conversion(seq2)

    mapped = []
    for idx, is_fwd in zip([idx_mk, idx_km], [True, False]):
        for hit1, hit2 in find_properly_paired_hits(
            filter_hits(
                idx.map(seq1_conv, seq2=seq2_conv, cs=False, MD=False), seq1, seq2
            ),
            fwd=is_fwd,
        ):
            tlen = max(hit1.r_en, hit2.r_en) - min(hit1.r_st, hit2.r_st)
            # flag1 and falg2 and both properly paired
            if is_fwd:
                flag1, flag2 = 99, 147
            else:
                flag1, flag2 = 83, 163
            ref1 = idx0.seq(hit1.ctg, hit1.r_st, hit1.r_en)
            ref2 = idx0.seq(hit2.ctg, hit2.r_st, hit2.r_en)
            if is_fwd:
                s1 = seq1
                s2 = mp.revcomp(seq2)
            else:
                s1 = mp.revcomp(seq1)
                s2 = seq2

            # fix soft clip
            # https://github.com/lh3/minimap2/issues/356
            cigar_str1 = hit1.cigar_str
            cigar1 = hit1.cigar
            if hit1.q_st > 0:
                cigar_str1 = f"{hit1.q_st}S" + cigar_str1
                cigar1 = [[hit1.q_st, 4]] + cigar1
            if hit1.q_en < len(s1):
                cigar_str1 = cigar_str1 + f"{len(s1) - hit1.q_en}S"
                cigar1 = cigar1 + [[len(s1) - hit1.q_en, 4]]
            cigar_str2 = hit2.cigar_str
            cigar2 = hit2.cigar
            if hit2.q_st > 0:
                cigar_str2 = f"{hit2.q_st}S" + cigar_str2
                cigar2 = [[hit2.q_st, 4]] + cigar2
            if hit2.q_en < len(s2):
                cigar_str2 = cigar_str2 + f"{len(s2) - hit2.q_en}S"
                cigar2 = cigar2 + [[len(s2) - hit2.q_en, 4]]

            md1, *debug = cal_md_and_tag(cigar1, s1, ref1, fwd_lib == is_fwd)
            md2, *_ = cal_md_and_tag(cigar2, s2, ref2, fwd_lib == is_fwd)
            # print(dict(zip(["Yf", "Zf", "Yc", "Zc", "NS", "NC"], debug)))
            # add this info into string and append ito name to debuging
            tag_info = f"__Yf:{debug[0]}+Zf:{debug[1]}+Yc:{debug[2]}+Zc:{debug[3]}+NS:{debug[4]}+NC:{debug[5]}"

            map1 = [
                name + tag_info,
                flag1,
                hit1.ctg,
                hit1.r_st + 1,
                hit1.mapq,
                cigar_str1,
                hit2.ctg,
                hit2.r_st + 1,
                tlen,
                s1,
                qua1,
                "MD:Z:" + md1,
                # strand tag
                f"ST:i:{int(is_fwd)}",
                # "cs:Z:" + hit1.cs,
            ]
            map2 = [
                name + tag_info,
                flag2,
                hit2.ctg,
                hit2.r_st + 1,
                hit2.mapq,
                cigar_str2,
                hit1.ctg,
                hit1.r_st + 1,
                -tlen,
                s2,
                qua2,
                "MD:Z:" + md2,
                f"ST:i:{int(is_fwd)}",
                # "cs:Z:" + hit2.cs,
            ]
            score = hit1.mapq + hit2.mapq
            mapped.append([score, map1, map2])

    random.shuffle(mapped)
    mapped = sorted(mapped, key=lambda x: x[0], reverse=True)
    for i, (_, map1, map2) in enumerate(mapped):
        if i > 0:
            map1[1] += 256
            map2[1] += 256
        print("\t".join(map(str, map1)))
        print("\t".join(map(str, map2)))


def map_file(ref_file, r1_file, r2_file, fwd_lib=True):
    # if mk and km file is not exist
    mk_file = ref_file + ".mk.fa"
    km_file = ref_file + ".km.fa"
    if not os.path.exists(mk_file) or not os.path.exists(km_file):
        mk_convert_file(ref_file, mk_file, km_file, include_ys_tag=False)

    # from .conversion import

    # s = a.seq("rRNA-Hsa-nucleus_locus", 100, 200)  # retrieve a subsequence from the index
    # print(mp.revcomp(s))  # reverse complement

    # Only map to the forward strand of the reference sequences. For paired-end
    # reads in the forward-reverse orientation, the first read is mapped to forward
    # strand of the reference and the second read to the reverse stand.
    # define MM_F_FOR_ONLY      (0x100000LL)
    # define MM_F_REV_ONLY      (0x200000LL)
    idx0 = mp.Aligner(
        fn_idx_in=ref_file,
        preset="sr",
        n_threads=8,
        k=10,
        w=10,
        min_cnt=0,
        min_chain_score=0,
        best_n=50,
    )
    idx_mk = mp.Aligner(
        fn_idx_in=mk_file,
        preset="sr",
        n_threads=8,
        k=10,
        w=10,
        min_cnt=0,
        min_chain_score=0,
        best_n=50,
        extra_flags=0x100000 if fwd_lib else 0x200000,
    )
    idx_km = mp.Aligner(
        fn_idx_in=km_file,
        preset="sr",
        n_threads=8,
        k=10,
        w=10,
        min_cnt=0,
        min_chain_score=0,
        best_n=50,
        extra_flags=0x200000 if fwd_lib else 0x100000,
    )
    # write sam header
    print("@HD\tVN:1.6\tSO:unsorted")
    # with qual
    for name, seq, *_ in mp.fastx_read(ref_file):
        print(f"@SQ\tSN:{name}\tLN:{len(seq)}")

    for (name1, seq1, qua1), (name2, seq2, qua2) in zip(
        mp.fastx_read(r1_file), mp.fastx_read(r2_file)
    ):
        if name1 != name2:
            raise ValueError("r1 and r2 not in the same order")
        run_mapping(name1, seq1, seq2, qua1, qua2, idx0, idx_km, idx_mk, fwd_lib)


if __name__ == "__main__":
    ref_file = "./test/ref.fa"
    r1_file = "./test/r1.fq.gz"
    r2_file = "./test/r2.fq.gz"
    map_file(ref_file, r1_file, r2_file, fwd_lib=False)

# # ## # on forward strand example
# seq_name = "test1"
# seq1 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACT"
# # seq1 = "TTGGGTTGTTTGGGGGTGTGGTGTGGGGTGGGT"
# qua1 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
# # seq2 = "AGTTTACCACCCGCTTTGGGCTGCATTCCCAAGCA"
# seq2 = "AGCCCACCACCCGCCCCGGGCCGCACCCCCAAGCA"
# qua2 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
#
# run_mapping(seq_name, seq1, seq2, qua1, qua2, idx0, idx_km, idx_mk, fwd_lib=True)
## # on rev strand example
## # seq1 = "GGGTCCTAACACGTGCGCTCGTGCTC"
## # seq2 = "TCTCGCCCGCCGCGCCGGGGAGGTGGAGCACGAG"
