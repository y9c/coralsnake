import cProfile
import pstats
from coralsnake.mapping import _map_batch_worker, _init_worker, _build_and_check_indices
import os
import multiprocessing as mp


def main():
    import pysam
    from bwamem import read_paired_fastx
    import time

    # Init worker
    ref_files = ["tests/benchmark/genes.fa", "tests/benchmark/transcript.fa"]
    i_dirs = [
        "tests/benchmark/test_idx_only/genes_idx",
        "tests/benchmark/test_idx_only/transcript_idx",
    ]
    ref_indices = [
        (None, None, str(i + 1) if len(ref_files) > 1 else "")
        for i in range(len(ref_files))
    ]
    ref_indices = _build_and_check_indices(ref_files, ref_indices, i_dirs, True)

    _init_worker(ref_indices, None, True, [{}, {}])  # dummy global rid maps

    it = read_paired_fastx(
        "tests/benchmark/trimmed_R1_100k.fq.gz", "tests/benchmark/trimmed_R2_100k.fq.gz"
    )
    batch = []
    for _ in range(2000):
        r1, r2 = next(it)
        batch.append(
            ((r1.name, r1.sequence, r1.quality), (r2.name, r2.sequence, r2.quality))
        )

    print("Profiling worker...")
    profiler = cProfile.Profile()
    profiler.enable()
    _map_batch_worker(batch, True, 6, 20, 0.8, 1.0, 1.0)
    profiler.disable()

    stats = pstats.Stats(profiler).sort_stats("tottime")
    stats.print_stats(20)


if __name__ == "__main__":
    main()
