import pysam

bam = pysam.AlignmentFile("tests/benchmark/verify_final_perf.bam", "rb")
ori1_primary = 0
ori2_secondary = 0
for read in bam:
    if not read.is_secondary:
        if read.get_tag("ST") == 1:
            ori1_primary += 1
    else:
        if read.get_tag("ST") == 2:
            ori2_secondary += 1
print(f"Ori1 primary: {ori1_primary}, Ori2 secondary: {ori2_secondary}")
