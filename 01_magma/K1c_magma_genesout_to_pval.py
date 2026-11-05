#!/usr/bin/env python3
"""Convert MAGMA .genes.out -> scDRS pval_file (.tsv: GENE + trait p-value, gene symbols).

MAGMA uses entrez IDs (from NCBI37.3.gene.loc); map entrez -> symbol.
"""
GENES_OUT = "/workspace/data/magma/finngen_O15.genes.out"
GENE_LOC = "/workspace/data/annotation/NCBI37.3.gene.loc"
OUT = "/workspace/data/magma/finngen_O15.pval.tsv"
TRAIT = "O15_PREECLAMPS"

# entrez -> symbol from hg19 gene loc
entrez2sym = {}
with open(GENE_LOC) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6:
            entrez2sym[f[0]] = f[-1]

rows = []
with open(GENES_OUT) as fh:
    header = None
    for line in fh:
        if line.startswith("#"):
            header = line.lstrip("#").split()
            continue
        if header is None:
            header = line.split()
            continue
        f = line.split()
        d = dict(zip(header, f))
        ent = d["GENE"]
        p = d["P"]
        sym = entrez2sym.get(ent)
        if not sym or sym == "":
            continue
        rows.append((sym, p))

with open(OUT, "w") as out:
    out.write(f"GENE\t{TRAIT}\n")
    for sym, p in rows:
        out.write(f"{sym}\t{p}\n")
print(f"[convert] genes mapped={len(rows)} -> {OUT}")
