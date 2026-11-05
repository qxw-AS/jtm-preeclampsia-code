#!/usr/bin/env python3
"""Integrate GSE173193 (memory-lean under a 4GB cgroup limit).

Run in two isolated processes to keep peak RSS low:

  python3 03_integrate_gse173193.py norm   # stage 1
      read raw counts -> normalize -> log1p -> HVG
      write combined_lognorm.h5ad (full genes) and combined_hvg.h5ad (HVG subset)

  python3 03_integrate_gse173193.py embed  # stage 2
      read combined_hvg.h5ad -> scale -> PCA -> Harmony(batch) -> neighbors
      -> leiden -> UMAP ; write combined_emb.h5ad

Raw counts remain in combined_raw.h5ad for scDRS.
"""
import gc
import sys

import scanpy as sc
import anndata as ad
import numpy as np

BASE = "/workspace/data/gse173193"
IN = f"{BASE}/combined_raw.h5ad"
LOGNORM = f"{BASE}/combined_lognorm.h5ad"
HVG_F = f"{BASE}/combined_hvg.h5ad"
EMB = f"{BASE}/combined_emb.h5ad"

sc.settings.verbosity = 1


def stage_norm():
    np.random.seed(0)
    a = sc.read_h5ad(IN)
    print(f"[load] {a.n_obs} cells x {a.n_vars} genes; batches={a.obs['batch'].nunique()}", flush=True)

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=3000, batch_key="batch", flavor="seurat")
    hvg = a.var.highly_variable.values
    print(f"[hvg] {int(hvg.sum())} HVGs", flush=True)

    a.write(LOGNORM)
    print(f"[write] {LOGNORM}", flush=True)

    ah = ad.AnnData(X=a.X[:, hvg], obs=a.obs.copy(), var=a.var.loc[hvg].copy())
    ah.write(HVG_F)
    print(f"[write] {HVG_F}  shape={ah.shape}", flush=True)


def stage_embed():
    np.random.seed(0)
    sc.settings.n_jobs = 1
    ah = sc.read_h5ad(HVG_F)
    print(f"[load hvg] {ah.n_obs} cells x {ah.n_vars} genes", flush=True)

    # scale densifies to float64 (~1.3GB); downcast to float32 to stay under the cgroup limit
    sc.pp.scale(ah, max_value=10)
    ah.X = ah.X.astype(np.float32)
    gc.collect()
    print("[scale] dense float32", flush=True)

    sc.tl.pca(ah, n_comps=30, svd_solver="arpack")
    sc.external.pp.harmony_integrate(ah, key="batch", max_iter_harmony=20)
    print("[harmony] X_pca_harmony present:", "X_pca_harmony" in ah.obsm, flush=True)

    sc.pp.neighbors(ah, n_neighbors=15, use_rep="X_pca_harmony")
    sc.tl.leiden(ah, resolution=0.8, flavor="igraph", n_iterations=2, directed=False, key_added="leiden")
    sc.tl.umap(ah)
    print(f"[leiden] {ah.obs['leiden'].nunique()} clusters", flush=True)

    emb = ad.AnnData(obs=ah.obs[["leiden"]].copy())
    emb.obsm["X_pca_harmony"] = ah.obsm["X_pca_harmony"]
    emb.obsm["X_umap"] = ah.obsm["X_umap"]
    emb.obsp["connectivities"] = ah.obsp["connectivities"]
    emb.obsp["distances"] = ah.obsp["distances"]
    emb.write(EMB)
    print(f"[done] embeddings -> {EMB}  clusters={ah.obs['leiden'].nunique()}", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "norm":
        stage_norm()
    elif mode == "embed":
        stage_embed()
    else:
        raise SystemExit("usage: norm|embed")
