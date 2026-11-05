#!/usr/bin/env python3
"""Marker-based cell-type annotation for placenta scRNA-seq.

Usage: python3 06_annotate.py <dataset: gse173193|gse282038>

Loads the log-normalized object, scores curated placental marker sets, assigns each
Leiden cluster the cell type with the highest mean score, attaches the embeddings,
and writes an annotated object (obs + embeddings).  Scoring uses the full-gene
log-normalized matrix so no `.raw` copy is required.
"""
import sys
import gc

import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd

MARKERS = {
    "Syncytiotrophoblast": ["CGA", "CGB", "CGB3", "PSG1", "PSG3", "PSG4", "CSH1", "CSH2",
                            "ERVFRD-1", "GCM1", "MUC15", "CYP19A1", "TFAP2A"],
    "Cytotrophoblast": ["KRT7", "GATA3", "TEAD4", "EGFR", "MSX2", "TP63", "PEG10", "HAND1"],
    "Extravillous_trophoblast": ["HLA-G", "ITGA1", "MMP2", "PAPPA2", "ASCL2", "NOTUM",
                                 "MMP11", "ADAM12", "PLAC8"],
    "Stromal_mesenchymal": ["COL1A1", "COL3A1", "COL1A2", "COL5A1", "DCN", "LUM",
                            "PDGFRA", "ACTA2", "FN1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "ENG", "PLVAP", "FLT1", "CD34"],
    "Macrophage_Hofbauer": ["CD68", "CD14", "CSF1R", "LYZ", "CD163", "MRC1", "MSR1", "AIF1"],
    "T_NK": ["CD3D", "CD3E", "NKG7", "KLRD1", "GNLY", "CCL5", "CD2", "TRAC"],
    "B_Plasma": ["CD79A", "MS4A1", "IGHG1", "MZB1", "JCHAIN", "CD79B"],
    "Erythrocyte_fetal": ["HBB", "HBA1", "HBA2", "ALAS2"],
}


def annotate(adata):
    """adata must carry full-gene log-normalized X and a 'leiden' column."""
    for name, genes in MARKERS.items():
        genes = [g for g in genes if g in adata.var_names]
        sc.tl.score_genes(adata, genes, score_name=f"score_{name}", ctrl_size=50)
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    cl = adata.obs["leiden"].astype(str)
    df = adata.obs[score_cols].copy()
    df["cluster"] = cl.values
    mean_by_cluster = df.groupby("cluster").mean()
    assign = mean_by_cluster.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["cell_type"] = cl.map(assign).astype("category")
    print("[annotate] cluster -> cell_type:")
    for c, t in assign.items():
        n = int((cl == c).sum())
        print(f"  cluster {c:>3}: {t:32s} n={n}")
    return adata, mean_by_cluster


def run_gse173193():
    base = "/workspace/data/gse173193"
    a = sc.read_h5ad(f"{base}/combined_lognorm.h5ad")  # full-gene log-normalized
    emb = sc.read_h5ad(f"{base}/combined_emb.h5ad")
    a.obs["leiden"] = emb.obs["leiden"].reindex(a.obs_names).values
    a.obsm["X_pca_harmony"] = emb.obsm["X_pca_harmony"]
    a.obsm["X_umap"] = emb.obsm["X_umap"]
    a.obsp["connectivities"] = emb.obsp["connectivities"]
    del emb
    gc.collect()

    a.obs["sample_condition"] = a.obs["condition"].astype(str)
    a, m = annotate(a)

    out = ad.AnnData(obs=a.obs.copy())
    out.obsm["X_umap"] = a.obsm["X_umap"]
    out.obsm["X_pca_harmony"] = a.obsm["X_pca_harmony"]
    out.obsp["connectivities"] = a.obsp["connectivities"]
    out.write(f"{base}/combined_annotated.h5ad")
    m.to_csv(f"{base}/combined_annotated_cluster_scores.csv")
    print(f"[done] -> {base}/combined_annotated.h5ad  types={a.obs['cell_type'].nunique()}")


def run_gse282038():
    base = "/workspace/data/gse282038"
    a = sc.read_h5ad(f"{base}/adata_integrated.h5ad")
    ann = a.raw.to_adata() if a.raw is not None else a.copy()
    ann.obs["leiden"] = a.obs["leiden"].values
    ann, m = annotate(ann)
    a.obs["cell_type"] = ann.obs["cell_type"].values
    out = ad.AnnData(obs=a.obs.copy())
    for k in ("X_umap", "X_pca"):
        if k in a.obsm:
            out.obsm[k] = a.obsm[k]
    out.write(f"{base}/adata_annotated.h5ad")
    m.to_csv(f"{base}/adata_annotated_cluster_scores.csv")
    print(f"[done] -> {base}/adata_annotated.h5ad")


if __name__ == "__main__":
    ds = sys.argv[1]
    if ds == "gse173193":
        run_gse173193()
    elif ds == "gse282038":
        run_gse282038()
    else:
        raise SystemExit("unknown dataset")
