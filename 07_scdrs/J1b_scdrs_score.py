#!/usr/bin/env python3
"""Run scDRS scoring for GSE173193 while respecting the 4GB cgroup limit.

scdrs.util.load_h5ad internally normalizes + filters with several full-size
copies and peaks above the limit.  Here we load the already log-normalized
object directly and drive scdrs.preprocess / scdrs.score_cell ourselves.
"""
import gc

import scanpy as sc
import scdrs
import numpy as np
import pandas as pd

BASE = "/workspace/data/gse173193"
GS = "/workspace/data/magma/O15_PREECLAMPS.gs"
OUT = "/workspace/data/scdrs/gse173193"
TRAIT = "O15_PREECLAMPS"
N_CTRL = 1000

a = sc.read_h5ad(f"{BASE}/combined_lognorm.h5ad")          # log-normalized, full genes
ann = sc.read_h5ad(f"{BASE}/combined_annotated.h5ad")      # obs only
ann = ann[a.obs_names]
a.obs["cell_type"] = ann.obs["cell_type"].values
del ann
gc.collect()
print(f"[load] {a.shape}", flush=True)

# Technical covariates -> forces scDRS implicit covariate-correction mode (sparse,
# chunked), which keeps peak memory under the 4GB cgroup limit.  Also regresses
# out library-size / mitochondrial content, which is standard scDRS practice.
cov = pd.DataFrame(index=a.obs_names)
cov["const"] = 1.0
cov["log_n_genes"] = np.log1p(a.obs["n_genes_by_counts"].astype(float).values)
cov["pct_mt"] = a.obs["pct_counts_mt"].astype(float).values
print(f"[cov] {cov.shape}", flush=True)

scdrs.preprocess(a, cov=cov, adj_prop="cell_type", n_mean_bin=20, n_var_bin=20, copy=False)
print("[preprocess] done (implicit cov mode)", flush=True)

dict_gs = scdrs.util.load_gs(GS, src_species="hsapiens", dst_species="hsapiens",
                             to_intersect=a.var_names)
gene_list, gene_weights = dict_gs[TRAIT]
print(f"[gs] {TRAIT}: n_gene={len(gene_list)}", flush=True)

df = scdrs.score_cell(
    a, gene_list, gene_weight=gene_weights, ctrl_match_key="mean_var",
    n_ctrl=N_CTRL, weight_opt="vs", return_ctrl_norm_score=True, verbose=True,
)
import os
os.makedirs(OUT, exist_ok=True)
df.iloc[:, 0:6].to_csv(f"{OUT}/{TRAIT}.score.gz", sep="\t", index=True, compression="gzip")
df.to_csv(f"{OUT}/{TRAIT}.full_score.gz", sep="\t", index=True, compression="gzip")
v_fdr = scdrs.util.multipletests(df["pval"].values, method="fdr_bh")[1]
print(f"[done] {TRAIT}: {int((v_fdr<0.1).sum())} FDR<0.1, {int((v_fdr<0.2).sum())} FDR<0.2 / {len(df)} cells", flush=True)
