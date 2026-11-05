#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Overlap correction of the replication cohort (checklist item F2).

The replication GWAS (GCST90269903) contains a FinnGen component, so the two
cohort estimates share correlated sampling error.  As a sensitivity analysis we
inflate replication-cohort standard errors by sqrt(296824 / 209254) = 1.19,
hold the between-cohort variance tau^2 at its primary estimate, and recompute
the DerSimonian-Laird random-effects meta-analysis.

Input
-----
results/figure4/fig4_mr_gene.csv   (output of 03_mr_meta_coloc)
    required columns: symbol, b_fin, se_fin, b_g2, se_g2, is_candidate, meta_ok

Output
------
results/figure4/fig4_overlap_correction.csv
    one row per gene with the primary and the inflated meta estimate
stdout
    counts of candidates with meta-P < 0.05 before and after inflation

Reference values (316 candidates)
---------------------------------
primary meta-P < 0.05 : 64
inflated meta-P < 0.05: 59
genes losing significance: NCK2, APOBEC2, PM20D1, HLA-DRB5, CPNE1
none of the nine direction-consistent genes changes status
(largest inflated meta-P: ELP3, 0.037)
"""
import math
import os

import numpy as np
import pandas as pd

IN = os.environ.get("PE_MR_GENE", "/workspace/results/figure4/fig4_mr_gene.csv")
OUT = os.environ.get("PE_OVERLAP_OUT", "/workspace/results/figure4/fig4_overlap_correction.csv")

N_REPL, N_NONFIN = 296824, 209254
INFLATION = math.sqrt(N_REPL / N_NONFIN)          # 1.19


def norm_p(b, se):
    """two-sided normal P value"""
    if b is None or se is None or se <= 0 or np.isnan(b) or np.isnan(se):
        return np.nan
    return math.erfc(abs(b / se) / math.sqrt(2.0))


def primary_tau2(b1, s1, b2, s2):
    """tau^2 of the primary (uninflated) two-cohort meta-analysis.

    k = 1 -> 0; a single-cohort gene carries no between-cohort variance.
    """
    bs = [(b, s) for b, s in ((b1, s1), (b2, s2))
          if b is not None and s is not None and not np.isnan(b) and not np.isnan(s) and s > 0]
    if len(bs) < 2:
        return 0.0
    w = np.array([1.0 / s ** 2 for _, s in bs])
    b = np.array([b for b, _ in bs])
    b_fe = float(np.sum(w * b) / np.sum(w))
    Q = float(np.sum(w * (b - b_fe) ** 2))
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    return max(0.0, (Q - (len(bs) - 1)) / C) if C > 0 else 0.0


def meta(b1, s1, b2, s2, tau2, inflate_replication=False):
    """DL random-effects meta with tau^2 held fixed.

    The inflation is applied to the replication-cohort SE whether or not that
    cohort is the only one available for the gene (replication-only genes are
    affected, FinnGen-only genes are not).
    """
    if s2 is not None and not np.isnan(s2) and inflate_replication:
        s2 = s2 * INFLATION
    bs = [(b, s) for b, s in ((b1, s1), (b2, s2))
          if b is not None and s is not None and not np.isnan(b) and not np.isnan(s) and s > 0]
    if not bs:
        return np.nan, np.nan, np.nan
    if len(bs) == 1:
        b, s = bs[0]
        return b, s, norm_p(b, s)
    w = np.array([1.0 / (s ** 2 + tau2) for _, s in bs])
    b = np.array([b for b, _ in bs])
    est = float(np.sum(w * b) / np.sum(w))
    se = math.sqrt(1.0 / np.sum(w))
    return est, se, norm_p(est, se)


def main():
    df = pd.read_csv(IN)
    df = df[df["is_candidate"].astype(str).str.lower().isin(["true", "1"])].copy()

    rows = []
    for _, r in df.iterrows():
        b1, s1 = r.get("b_fin"), r.get("se_fin")
        b2, s2 = r.get("b_g2"), r.get("se_g2")
        tau2 = primary_tau2(b1, s1, b2, s2)
        bp, sp, pp = meta(b1, s1, b2, s2, tau2, inflate_replication=False)
        ba, sa, pa = meta(b1, s1, b2, s2, tau2, inflate_replication=True)
        rows.append(dict(gene=r["symbol"], ensg=r.get("ensg"), tau2=tau2,
                         b_meta_primary=bp, se_meta_primary=sp, p_meta_primary=pp,
                         b_meta_inflated=ba, se_meta_inflated=sa, p_meta_inflated=pa,
                         sig_primary=(pp < 0.05), sig_inflated=(pa < 0.05)))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)

    n = len(out)
    n0 = int(out["sig_primary"].sum())
    n1 = int(out["sig_inflated"].sum())
    lost = sorted(out.loc[out["sig_primary"] & ~out["sig_inflated"], "gene"])

    print(f"[overlap] candidates = {n}")
    print(f"[overlap] inflation factor = {INFLATION:.3f}")
    print(f"[overlap] meta-P < 0.05: {n0} -> {n1}")
    print(f"[overlap] losing significance ({len(lost)}): {', '.join(lost)}")
    print(f"[overlap] written -> {OUT}")


if __name__ == "__main__":
    main()
