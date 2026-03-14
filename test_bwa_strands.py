from bwamem import BwaAligner

a = BwaAligner("tests/benchmark/test_idx_only/genes_idx/ref.mk", min_seed_len=14, mark_secondary=True, softclip_supplementary=True)
print(a.align_raw("GTTGGTTGGTGGTGTTGTTGGTTGGTTGGTGGTGTTGTTG"))
print(a.align_raw("CAACCAACCACCACAACAACCAACCAACCACCACAACAAC"))
