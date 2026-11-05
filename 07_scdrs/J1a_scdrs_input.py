#!/usr/bin/env python3
"""Build the scDRS input for GSE173193: raw counts + cell-type / condition metadata.

scDRS `compute_score` expects raw counts in .X (default --flag_raw_count=True applies
its own size-factor normalization + log1p).  We attach the cell_type / condition from
the annotated object (matched by obs_names) so `--adj_prop cell_type` can be used.
"""
import gc

import scanpy as sc
import anndata as ad

BASE = "/workspace/data/gse173193"
a = sc.read_h5ad(f"{BASE}/combined_raw.h5ad")          # raw counts
ann = sc.read_h5ad(f"{BASE}/combined_annotated.h5ad")  # obs only

ann = ann[a.obs_names]  # align
a.obs["cell_type"] = ann.obs["cell_type"].values
a.obs["sample_condition"] = a.obs["condition"].astype(str)

out = f"{BASE}/combined_scdrs.h5ad"
a.write(out)
print(f"[done] {out}  shape={a.shape}  types={a.obs['cell_type'].nunique()}")

gc.collect()
