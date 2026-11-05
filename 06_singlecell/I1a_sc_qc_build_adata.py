#!/usr/bin/env python3
"""Build per-sample AnnData for GSE173193 (placenta, 10x) + GSE282038 (EOPE trophoblast),
run per-sample QC, and concatenate into one object per dataset.

Usage: python3 01_build_adata.py <dataset: gse173193|gse282038>
"""
import os
import sys
import glob
import shutil
import tempfile

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

MT_PREFIX = "MT-"


def read_one(sample_dir_files, sample_id, condition):
    """Given a dict {barcodes,features,matrix} paths, build AnnData with proper 10x naming."""
    tmp = tempfile.mkdtemp(prefix=f"10x_{sample_id}_")
    try:
        shutil.copy(sample_dir_files["barcodes"], os.path.join(tmp, "barcodes.tsv.gz"))
        shutil.copy(sample_dir_files["features"], os.path.join(tmp, "features.tsv.gz"))
        shutil.copy(sample_dir_files["matrix"], os.path.join(tmp, "matrix.mtx.gz"))
        a = sc.read_10x_mtx(tmp, var_names="gene_symbols", cache=False)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    a.var_names_make_unique()
    a.obs["sample"] = sample_id
    a.obs["condition"] = condition
    return a


def qc(a, min_genes=200, max_pct_mt=20.0, min_cells=3):
    a.var["mt"] = a.var_names.str.startswith(MT_PREFIX)
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    n0 = a.n_obs
    a = a[a.obs["n_genes_by_counts"] >= min_genes].copy()
    a = a[a.obs["pct_counts_mt"] <= max_pct_mt].copy()
    sc.pp.filter_genes(a, min_cells=min_cells)
    print(f"  QC {a.obs['sample'].iloc[0] if a.n_obs else '?'}: {n0} -> {a.n_obs} cells, {a.n_vars} genes")
    return a


def build_gse173193():
    base = "/workspace/data/gse173193"
    rows = []
    with open(os.path.join(base, "filelist.txt")) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            gsm, label, group, url = line.split("\t")
            fn = url.split("/")[-1]
            rows.append((gsm, label, group, fn))
    df = pd.DataFrame(rows, columns=["gsm", "label", "group", "fn"])

    # Write each sample to disk first (memory-light), then concatenate on disk.
    sample_dir = os.path.join(base, "samples")
    os.makedirs(sample_dir, exist_ok=True)
    per_sample = []
    labels = []
    for (gsm, label, group), sub in df.groupby(["gsm", "label", "group"]):
        files = {}
        for _, r in sub.iterrows():
            if r["fn"].endswith("barcodes.tsv.gz"):
                files["barcodes"] = os.path.join(base, r["fn"])
            elif r["fn"].endswith("features.tsv.gz"):
                files["features"] = os.path.join(base, r["fn"])
            elif r["fn"].endswith("matrix.mtx.gz"):
                files["matrix"] = os.path.join(base, r["fn"])
        if len(files) != 3:
            print(f"  skip {label}: incomplete files {list(files)}")
            continue
        print(f"[read] {label} ({group})")
        a = read_one(files, label, group)
        a = qc(a)
        a.obs["tsample"] = label
        a.obs["dataset"] = "GSE173193"
        a.X = a.X.tocsr()
        fp = os.path.join(sample_dir, f"{label}.h5ad")
        a.write(fp)
        per_sample.append(fp)
        labels.append(label)
        del a
    if not per_sample:
        raise SystemExit("no samples loaded for GSE173193")
    out = "/workspace/data/gse173193/combined_raw.h5ad"
    ad.experimental.concat_on_disk(
        per_sample, out, join="outer", label="batch", keys=labels, index_unique="-"
    )
    print(f"[done] GSE173193 combined -> {out}")


def build_gse282038():
    base = "/workspace/data/gse282038"
    gsm = "GSM8634701"
    files = {
        "barcodes": os.path.join(base, f"{gsm}_barcodes.tsv.gz"),
        "features": os.path.join(base, f"{gsm}_features.tsv.gz"),
        "matrix": os.path.join(base, f"{gsm}_matrix.mtx.gz"),
    }
    for k, v in files.items():
        if not os.path.exists(v):
            raise SystemExit(f"missing {v}")
    a = read_one(files, "EOPE_FGR", "EOPE")
    a = qc(a)
    a.obs["dataset"] = "GSE282038"
    out = "/workspace/data/gse282038/adata_raw.h5ad"
    a.write(out)
    print(f"[done] GSE282038: {a.n_obs} cells x {a.n_vars} genes -> {out}")
    return a


if __name__ == "__main__":
    ds = sys.argv[1]
    if ds == "gse173193":
        build_gse173193()
    elif ds == "gse282038":
        build_gse282038()
    else:
        raise SystemExit("unknown dataset")
