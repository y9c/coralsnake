import multiprocessing as mp
from bwamem import BwaAligner
import time

opts = {"min_seed_len": 14, "mark_secondary": True, "softclip_supplementary": True}

print("Loading index in main process...")
t0 = time.time()
global_a = BwaAligner("tests/benchmark/test_idx_only/genes_idx/ref.mk", **opts)
print(f"Loaded in {time.time() - t0:.2f}s")


def worker_func(seq):
    # Use the global_a which was inherited via fork
    res = global_a.align_raw(seq)
    return len(res)


if __name__ == "__main__":
    t0 = time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(4) as pool:
        seqs = ["ACGTACGTACGTACGT" * 5] * 10000
        pool.map(worker_func, seqs)
    print(f"Workers done in {time.time() - t0:.2f}s")
