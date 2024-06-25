#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-08 20:32


import os

import mappy as mp

# on forward strand example
seq1 = "TCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGT"
seq2 = "AGTTTACCACCCGCTTTGGGCTGCATTCCCAAGCA"

idx = mp.Aligner(fn_idx_in="./debug.fa", preset="sr")

# q_st  q_en  strand  ctg  ctg_len  r_st  r_en  mlen  blen  mapq  cg:Z:cigar_str

for hit in idx.map(seq1, seq2=seq2, cs=True, MD=True):
    print(str(hit))
print(".......")

# reverse
for hit in idx.map(seq1, seq2=seq2[::-1], cs=True, MD=True):
    print(str(hit))
print(".......")

# complement
for hit in idx.map(seq1, seq2=mp.revcomp(seq2)[::-1], cs=True, MD=True):
    print(str(hit))
print(".......")

# rev complement complement
for hit in idx.map(seq1, seq2=mp.revcomp(seq2), cs=True, MD=True):
    print(str(hit))
print(".......")
