#!/usr/bin/env python3
"""scDRS downstream analyses: group-level cell-type association (Runs A/B/C),
donor-level aggregation, for GSE173193 and GSE282038.

Run A = all cells; Run B = control-donor cells only; Run C = PE-donor cells only.
MC association test = 95th percentile of normalized disease score vs null control scores.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import scanpy as sc
import scdrs


def mc_assoc(df_reg, cells, n_ctrl):
    """Group-level MC association: q95 of disease score vs control q95 distribution."""
    sub = df_reg.loc[cells]
    score_q95 = np.quantile(sub["norm_score"].values, 0.95)
    ctrl = sub[[c for c in df_reg.columns if c.startswith("ctrl_norm_score")]].values
    v_ctrl_q95 = np.quantile(ctrl, 0.95, axis=0)
    mc_p = ((v_ctrl_q95 >= score_q95).sum() + 1) / (v_ctrl_q95.shape[0] + 1)
    mc_z = (score_q95 - v_ctrl_q95.mean()) / (v_ctrl_q95.std() + 1e-12)
    return mc_p, mc_z


def group_table(df_reg, group_col, adata):
    groups = sorted(adata.obs[group_col].astype(str).unique())
    rows = []
    for g in groups:
        cells = adata.obs_names[adata.obs[group_col].astype(str) == g]
        cells = [c for c in cells if c in df_reg.index]
        if len(cells) < 10:
            continue
        p_fdr = scdrs.util.multipletests(df_reg.loc[cells, "pval"].values, method="fdr_bh")[1]
        p, z = mc_assoc(df_reg, cells, None)
        rows.append({
            "group": g, "n_cell": len(cells),
            "n_fdr_0.05": int((p_fdr < 0.05).sum()),
            "n_fdr_0.1": int((p_fdr < 0.1).sum()),
            "n_fdr_0.2": int((p_fdr < 0.2).sum()),
            "prop_fdr_0.1": float((p_fdr < 0.1).mean()),
            "assoc_mcp": p, "assoc_mcz": z,
        })
    return pd.DataFrame(rows)


def run_dataset(adata_file, score_file, out_prefix, condition_key="condition"):
    adata = sc.read_h5ad(adata_file)
    # map full-score cells to adata cells
    dict_score = scdrs.util.load_scdrs_score(score_file, obs_names=adata.obs_names)
    trait = list(dict_score.keys())[0]
    df_full = dict_score[trait]
    # intersect
    common = sorted(set(adata.obs_names) & set(df_full.index))
    adata = adata[common].copy()
    df_full = df_full.loc[common]
    # normalize score to a df_reg with norm_score + ctrl columns
    df_reg = df_full.copy()

    res = {}
    # Run A: all cells
    res["runA_cell_type"] = group_table(df_reg, "cell_type", adata)
    if condition_key in adata.obs:
        res["runA_condition"] = group_table(df_reg, condition_key, adata)
    # Run B / C for GSE173193 (PE vs control donors)
    if "sample_condition" in adata.obs:
        for tag, cond in [("runB_control", "control group"), ("runC_PE", "PE group")]:
            sub = adata.obs_names[adata.obs["sample_condition"].astype(str) == cond]
            sub = [c for c in sub if c in df_reg.index]
            if len(sub) > 0:
                res[tag + "_cell_type"] = group_table(df_reg, "cell_type", adata[sub])
    # donor-level aggregation
    if "sample" in adata.obs:
        df_donor = df_reg.loc[common, ["norm_score"]].copy()
        df_donor["sample"] = adata.obs.loc[common, "sample"].values
        if "sample_condition" in adata.obs:
            df_donor["sample_condition"] = adata.obs.loc[common, "sample_condition"].values
        agg = df_donor.groupby("sample").agg(
            mean_norm_score=("norm_score", "mean"),
            median_norm_score=("norm_score", "median"),
            n_cell=("norm_score", "size"),
        )
        if "sample_condition" in df_donor:
            agg["condition"] = df_donor.groupby("sample")["sample_condition"].first()
        res["donor"] = agg.reset_index()

    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    for k, v in res.items():
        v.to_csv(f"{out_prefix}.{k}.csv", index=False)
        print(f"[out] {out_prefix}.{k}.csv  ({v.shape[0]} rows)")
    return res


if __name__ == "__main__":
    ds = sys.argv[1]
    if ds == "gse173193":
        run_dataset("/workspace/data/gse173193/combined_annotated.h5ad",
                    "/workspace/data/scdrs/gse173193/O15_PREECLAMPS.full_score.gz",
                    "/workspace/results/gse173193")
    elif ds == "gse282038":
        run_dataset("/workspace/data/gse282038/adata_annotated.h5ad",
                    "/workspace/data/scdrs/gse282038/O15_PREECLAMPS.full_score.gz",
                    "/workspace/results/gse282038")
