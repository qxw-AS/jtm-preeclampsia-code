#!/usr/bin/env python3
"""Build MAGMA --pval input from FinnGen R12 O15_PREECLAMPS summary stats.

Output columns: SNP, P, N  (MAGMA: --pval finngen.pval use=SNP,P ncol=N)
FinnGen 'rsids' may contain comma-separated rs IDs; take the first; skip if none.
"""
import gzip

SRC = "/workspace/data/finngen/finngen_R12_O15_PREECLAMPS.gz"
DST = "/workspace/data/finngen/finngen_O15_PREECLAMPS.magma.pval"
N = 268336  # 9023 cases + 259313 controls

n_in = 0
n_out = 0
with gzip.open(SRC, "rt") as fh, open(DST, "w") as out:
    header = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    i_rs = idx["rsids"]
    i_p = idx["pval"]
    out.write("SNP\tP\tN\n")
    for line in fh:
        n_in += 1
        f = line.rstrip("\n").split("\t")
        rs = f[i_rs]
        if (not rs) or rs == "NA":
            continue
        rs = rs.split(",")[0]
        p = f[i_p]
        if not p or p == "NA":
            continue
        out.write(f"{rs}\t{p}\t{N}\n")
        n_out += 1
print(f"[pval] read={n_in} written={n_out} -> {DST}")
