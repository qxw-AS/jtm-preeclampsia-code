#!/usr/bin/env bash
# MAGMA gene-level analysis of FinnGen O15_PREECLAMPS using the 1000G EUR reference (g1000_eur, hg19)
# and hg19 gene locations (NCBI37.3.gene.loc). 10-kb window SNP->gene mapping (scDRS paper setting).
set -euo pipefail
MAGMA=/workspace/work/magma/magma
REF=/workspace/data/ref1000g/g1000_eur
GENE_LOC=/workspace/data/annotation/NCBI37.3.gene.loc
PVAL=/workspace/data/finngen/finngen_O15_PREECLAMPS.magma.pval
OUT=/workspace/data/magma
mkdir -p "$OUT"

echo "== [1/2] annotate (window=10,10) =="
"$MAGMA" --annotate window=10,10 \
  --snp-loc "${REF}.bim" \
  --gene-loc "$GENE_LOC" \
  --out "$OUT/g1000_eur"

echo "== [2/2] gene analysis =="
"$MAGMA" --bfile "$REF" \
  --pval "$PVAL" use=SNP,P ncol=N \
  --gene-annot "$OUT/g1000_eur.genes.annot" \
  --out "$OUT/finngen_O15"

echo "== done: $OUT/finngen_O15.genes.out =="
wc -l "$OUT/finngen_O15.genes.out"
