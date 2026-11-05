#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Figure 5 — 数据准备 (Layer 2 溯源, 全部真实数据重算)
=====================================================================================
MR 全套敏感性诊断 (V2.3 Suppl Fig5): Cochran's Q、MR-Egger 截距、
漏斗图、留一法分析、Steiger 检验汇总。

数据源 (MR 流水线中间产物):
  · 工具变量(暴露): results/figure4/fig4_instruments.csv
        Peng_2017 胎盘 cis-eQTL SuSiE 可信集 (n=149); 每 IV: b_exp,se_exp,F,rsid
  · 结局-1: results/figure4/fig4_outcome_finngen.csv (FinnGen R12 O15_PREECLAMPS)
        b_out_h 已对齐到暴露效应等位; N_case=9023, N_control=259313
  · 样本量: N_eqtl=149, N_fin=268336
方法 (透明、可溯源):
  · 单 IV MR:  Wald ratio = b_out/b_exp, se = se_out/|b_exp|
  · ≥2 IV MR:  IVW(固定效应), Cochran's Q, MR-Egger(WLS b_out~b_exp, 权重 1/se_out^2)
  · 留一法:    逐个剔除 IV 重算 IVW
  · Steiger(近似):  r2 = z^2/(z^2+N); 正确方向 = r2_exposure > r2_outcome
输出 -> /workspace/results/suppl5/
"""
import json, os, math
import numpy as np
import pandas as pd
from scipy import stats

BASE = "/workspace"
OUT = f"{BASE}/results/suppl5"; os.makedirs(OUT, exist_ok=True)
F4 = f"{BASE}/results/figure4"

N_EQTL = 149
N_FIN = 268336

inst = pd.read_csv(f"{F4}/fig4_instruments.csv")
inst["rsid"] = inst["rsid"].astype(str)
fin = pd.read_csv(f"{F4}/fig4_outcome_finngen.csv")
fin["rsid"] = fin["rsid"].astype(str)
fin_map = (fin.sort_values("p_out").drop_duplicates("rsid")
           .set_index("rsid")[["b_out_h", "se_out", "b_out", "af_alt"]])
df = inst.merge(fin_map, left_on="rsid", right_index=True, how="left")
df = df[df["b_out_h"].notna()].copy()
print(f"[merge] IVs with FinnGen outcome: {len(df)} across {df['symbol'].nunique()} genes")


def z2r2(z, N):
    return z * z / (z * z + N)


def wls_egger(b_exp, b_out, se_out):
    """Weighted regression b_out = a + slope*b_exp, weights=1/se_out^2. Return a,se_a,slope,p_a."""
    w = 1.0 / se_out ** 2
    Sw = w.sum(); Swx = (w * b_exp).sum(); Swy = (w * b_out).sum()
    Swxx = (w * b_exp * b_exp).sum(); Swxy = (w * b_exp * b_out).sum()
    denom = Sw * Swxx - Swx ** 2
    if denom <= 0:
        return np.nan, np.nan, np.nan, np.nan
    slope = (Sw * Swxy - Swx * Swy) / denom
    a = (Swy - slope * Swx) / Sw
    # intercept SE
    yhat = a + slope * b_exp
    dof = len(b_exp) - 2
    if dof <= 0:
        return a, np.nan, slope, np.nan
    resid = b_out - yhat
    s2 = (w * resid ** 2).sum() / dof
    se_a = math.sqrt(s2 * Swxx / denom) if denom > 0 else np.nan
    p_a = 2 * stats.t.sf(abs(a / se_a), dof) if se_a and not np.isnan(se_a) else np.nan
    return a, se_a, slope, p_a


rows, per_iv_rows = [], []
for ensg, g in df.groupby("ensg"):
    sym = g["symbol"].iloc[0]
    ne = g["n_iv"].iloc[0] if "n_iv" in g else len(g)
    b_exp = g["b_exp"].values; se_exp = g["se_exp"].values
    b_out = g["b_out_h"].values; se_out = g["se_out"].values
    ratio = b_out / b_exp; se_ratio = se_out / np.abs(b_exp)
    k = len(g)
    # IVW
    w = 1.0 / se_ratio ** 2
    b_ivw = np.sum(w * ratio) / np.sum(w); se_ivw = math.sqrt(1.0 / np.sum(w))
    p_ivw = 2 * stats.norm.sf(abs(b_ivw / se_ivw))
    Q = np.sum(w * (ratio - b_ivw) ** 2); dfq = k - 1
    Qp = stats.chi2.sf(Q, dfq) if dfq >= 1 else np.nan
    # Egger (needs >=2 IV)
    a, se_a, slope, p_a = (wls_egger(b_exp, b_out, se_out) if k >= 2 else (np.nan,) * 4)
    # LOO
    loo_flip = 0
    if k >= 2:
        for i in range(k):
            m = np.ones(k, bool); m[i] = False
            wi = 1.0 / se_ratio[m] ** 2
            bi = np.sum(wi * ratio[m]) / np.sum(wi)
            if np.sign(bi) != np.sign(b_ivw):
                loo_flip += 1
    rows.append(dict(ensg=ensg, symbol=sym, n_iv=k, b_ivw=b_ivw, se_ivw=se_ivw, p_ivw=p_ivw,
                     Q=Q, Q_df=dfq, Q_p=Qp, egger_intercept=a, egger_se=se_a, egger_p=p_a,
                     egger_slope=slope, n_loo_flip=loo_flip,
                     Fmin=g["F"].min(), Fmax=g["F"].max()))
    # per-IV (for funnel / LOO plots on representatives)
    for _, r in g.iterrows():
        z_exp = r["b_exp"] / r["se_exp"]
        z_out = r["b_out_h"] / r["se_out"]
        per_iv_rows.append(dict(symbol=sym, ensg=ensg, rsid=r["rsid"], b_exp=r["b_exp"],
                                se_exp=r["se_exp"], b_out=r["b_out_h"], se_out=r["se_out"],
                                ratio=r["b_out_h"] / r["b_exp"], se_ratio=r["se_out"] / abs(r["b_exp"]),
                                F=r["F"], pip_exp=r.get("pip_exp", np.nan),
                                r2_exp=z2r2(z_exp, N_EQTL), r2_out=z2r2(z_out, N_FIN)))

diag = pd.DataFrame(rows)
periv = pd.DataFrame(per_iv_rows)
periv["steiger_correct"] = periv["r2_exp"] > periv["r2_out"]

# ---- candidate / Tier 标注 (直接继承 fig4) ----
mr = pd.read_csv(f"{F4}/fig4_mr_gene.csv")[["ensg", "Tier", "PPH4", "is_candidate", "fdr_fin"]]
diag = diag.merge(mr, on="ensg", how="left")
periv = periv.merge(mr[["ensg", "Tier", "is_candidate"]], on="ensg", how="left")

diag.to_csv(f"{OUT}/suppl5_multi_iv_diag.csv", index=False)
periv.to_csv(f"{OUT}/suppl5_per_iv.csv", index=False)

# ---- representative genes for funnel / LOO ----
multi = diag.sort_values(["n_iv", "p_ivw"], ascending=[False, True])
reps = list(multi["symbol"].head(3))
meta = {
    "provenance": {
        "iv": "fig4_instruments.csv (Peng_2017 placenta cis-eQTL SuSiE, n=149)",
        "outcome": "fig4_outcome_finngen.csv (FinnGen R12 O15_PREECLAMPS, 9023/259313)",
        "methods": "IVW / Cochran Q / MR-Egger WLS / LOO / Steiger-approx r2=z2/(z2+N)",
    },
    "n_genes_total_iv": int(df["ensg"].nunique()),
    "n_multi_iv_genes": int(len(diag)),
    "n_iv_dist": diag["n_iv"].value_counts().sort_index().to_dict(),
    "n_egger_tested": int(diag["egger_p"].notna().sum()),
    "n_Q_signif_005": int((diag["Q_p"] < 0.05).sum()),
    "n_egger_intercept_signif_005": int((diag["egger_p"] < 0.05).sum()),
    "n_loo_flip_genes": int((diag["n_loo_flip"] > 0).sum()),
    "fstat_min": float(df["F"].min()), "fstat_median": float(df["F"].median()),
    "fstat_max": float(df["F"].max()),
    "n_iv_over_f10": int((df["F"] > 10).sum()), "n_iv_total": int(len(df)),
    "steiger_correct_frac": float(periv["steiger_correct"].mean()),
    "steiger_n": int(len(periv)),
    "reps": reps,
}
with open(f"{OUT}/suppl5_meta.json", "w") as fh:
    json.dump(meta, fh, indent=2, ensure_ascii=False)
print(json.dumps(meta, indent=2, ensure_ascii=False))
print("\n[reps]", reps)
print(diag[diag["symbol"].isin(reps)][["symbol", "n_iv", "b_ivw", "p_ivw", "Q_p", "egger_p", "Tier"]].to_string(index=False))
