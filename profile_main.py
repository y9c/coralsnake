import cProfile
import pstats
from coralsnake.cli import map

def main():
    print("Profiling main process...")
    profiler = cProfile.Profile()
    profiler.enable()
    
    try:
        map.callback(
            r1_file="tests/benchmark/trimmed_R1_100k.fq.gz",
            r2_file="tests/benchmark/trimmed_R2_100k.fq.gz",
            ref_files=["tests/benchmark/genes.fa", "tests/benchmark/transcript.fa"],
            output_files=["tests/benchmark/verify_final_perf.bam"],
            unmap_file=None,
            max_mismatches=6,
            threads=12,
            min_alignment_length=20,
            min_mapping_ratio=0.8,
            max_a2g_ratio=1.0,
            max_c2t_ratio=1.0,
            index_dir=["tests/benchmark/test_idx_only/genes_idx", "tests/benchmark/test_idx_only/transcript_idx"],
            index_only=False,
            batch_size=2000,
            library_type="forward",
            reference_strand="double"
        )
    except Exception as e:
        print(e)
    
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats('tottime')
    stats.print_stats(30)

if __name__ == "__main__":
    main()
