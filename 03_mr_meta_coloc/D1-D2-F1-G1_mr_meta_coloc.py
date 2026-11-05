#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure4 上游 (3/6): 两阶段三级 MR → 双 GWAS 随机效应 meta → COLOC Tier 分级 → 漏斗
=====================================================================================
数据源 (Layer 2 溯源):
  · 工具变量: fig4_instruments.csv (Peng_2017 胎盘 eQTL SuSiE 可信集, n=149)
  · 结局-1 : fig4_outcome_finngen.csv (FinnGen R12 O15_PREECLAMPS, 9023/259313)
  · 结局-2 : fig4_outcome_gwas2.csv (GCST90269903 Tyrmi2023 PE-meta, n=296824)
  · COLOC   : QTD000394.credible_sets.tsv.gz (胎盘 eQTL 可信集 β/se)
              + Finngen 全量汇总统计 (区域 β/se)  -> 近似 coloc.abf (ABF, 独立近似)

方法:
  · IV 统计量: 单 IV -> Wald ratio (b_out/b_exp, se_out/|b_exp|)
              ≥2 IV -> IVW (固定效应, w=1/se_ratio^2)
  · meta    : 两套 GWAS MR 效应随机效应 (DerSimonian-Laird)
  · 稳健判定: 方向一致 & meta-p<0.05
  · COLOC   : Wakefield ABF (coloc.abf 独立近似), p1=p2=1e-4, p12=1e-5
              Tier1 PPH4>0.80 / Tier2 0.50<PPH4<=0.80 / else 剔除
输出:
  · fig4_mr_gene.csv    每基因: 两 GWAS MR + meta + COLOC + Tier
  · fig4_funnel.csv     三级漏斗剩余基因数 (含 Tier1/Tier2)
"""
import os, gzip, math
import numpy as np
import pandas as pd
from scipy import stats

OUT = "/workspace/results/figure4"
FIN_FULL = "/workspace/data/finngen/finngen_R12_O15_PREECLAMPS.gz"
CRED = "/workspace/data/eqtl_placenta/QTD000394.credible_sets.tsv.gz"

# 样本量 (Layer2 溯源)
N_FIN, N_FIN_CASE = 268336, 9023            # FinnGen R12 O15_PREECLAMPS (9023 cases / 259313 controls)
N_G2 = 296824                                # GCST90269903 (Tyrmi 2023 PE meta)
N_EQTL = 149                                 # Peng_2017 胎盘 eQTL

inst = pd.read_csv(f"{OUT}/fig4_instruments.csv")
inst["rsid"] = inst["rsid"].astype(str)

# ---------------- 1) 结局效应表 ----------------
fin = pd.read_csv(f"{OUT}/fig4_outcome_finngen.csv")
fin["rsid"] = fin["rsid"].astype(str)
g2 = pd.read_csv(f"{OUT}/fig4_outcome_gwas2.csv")
g2["rsid"] = g2["rsid"].astype(str)

# finngen: b_out_h 已对齐到暴露 ea ; gwas2: b_out 已对齐到暴露 ea (见 21 脚本)
fin_map = (fin.sort_values("p_out").drop_duplicates("rsid")
           .set_index("rsid")[["b_out_h", "se_out", "af_alt"]]
           .rename(columns={"b_out_h": "b_fin", "se_out": "se_fin", "af_alt": "af_fin"}))
g2_map = (g2.sort_values("p_out").drop_duplicates("rsid")
          .set_index("rsid")[["b_out", "se_out"]]
          .rename(columns={"b_out": "b_g2", "se_out": "se_g2"}))

df = inst.merge(fin_map, left_on="rsid", right_index=True, how="left")
df = df.merge(g2_map, left_on="rsid", right_index=True, how="left")
print(f"[merge] IVs={len(df)}  fin matched={df['b_fin'].notna().sum()}  "
      f"g2 matched={df['b_g2'].notna().sum()}", flush=True)

def _ivw(b_exp, se_exp, b_out, se_out):
    """返回 (beta, se, npairs, method)。单 IV=Wald; ≥2=IVW(固定)。"""
    m = ~(np.isnan(b_exp) | np.isnan(b_out) | np.isnan(se_out) | np.isnan(se_exp))
    b_exp, se_exp = np.array(b_exp)[m], np.array(se_exp)[m]
    b_out, se_out = np.array(b_out)[m], np.array(se_out)[m]
    if len(b_exp) == 0:
        return np.nan, np.nan, 0, ""
    ratio = b_out / b_exp
    se_ratio = se_out / np.abs(b_exp)
    if len(b_exp) == 1:
        return ratio[0], se_ratio[0], 1, "Wald"
    w = 1.0 / se_ratio ** 2
    b = np.sum(w * ratio) / np.sum(w)
    se = math.sqrt(1.0 / np.sum(w))
    return b, se, len(b_exp), "IVW"

rows = []
for ensg, g in df.groupby("ensg"):
    sym = g["symbol"].iloc[0]
    b_f, se_f, np_f, m_f = _ivw(g["b_exp"], g["se_exp"], g["b_fin"], g["se_fin"])
    b_2, se_2, np_2, m_2 = _ivw(g["b_exp"], g["se_exp"], g["b_g2"], g["se_g2"])
    rows.append(dict(ensg=ensg, symbol=sym, n_iv=len(g),
                     b_fin=b_f, se_fin=se_f, n_iv_fin=np_f, method_fin=m_f,
                     b_g2=b_2, se_g2=se_2, n_iv_g2=np_2, method_g2=m_2))
mr = pd.DataFrame(rows)

def _pval(b, se):
    return 2 * stats.norm.sf(np.abs(b / se))

mr["p_fin"] = _pval(mr["b_fin"].values, mr["se_fin"].values)
mr["p_g2"] = _pval(mr["b_g2"].values, mr["se_g2"].values)

# ---------------- 2) 双 GWAS 随机效应 meta (DerSimonian-Laird) ----------------
def dl_meta(b1, se1, b2, se2):
    if np.isnan(b1) and np.isnan(b2):
        return np.nan, np.nan, np.nan, False
    bs = np.array([b1, b2], float); ses = np.array([se1, se2], float)
    m = ~np.isnan(bs)
    bs, ses = bs[m], ses[m]
    if len(bs) == 1:
        # 单队列回退: 无第二队列时原样返回该队列估计, 便于下游展示,
        # 但这不是 meta 分析, 故 meta_ok = False (k = 1).
        b, se = bs[0], ses[0]
        return b, se, _pval(b, se), False
    w = 1.0 / ses ** 2
    b_fe = np.sum(w * bs) / np.sum(w)
    Q = np.sum(w * (bs - b_fe) ** 2)
    k = len(bs)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1.0 / (ses ** 2 + tau2)
    b = np.sum(wr * bs) / np.sum(wr)
    se = math.sqrt(1.0 / np.sum(wr))
    return b, se, _pval(b, se), True

mr[["b_meta", "se_meta", "p_meta", "meta_ok"]] = mr.apply(
    lambda r: pd.Series(dl_meta(r["b_fin"], r["se_fin"], r["b_g2"], r["se_g2"])), axis=1)
mr["dir_consistent"] = np.sign(mr["b_fin"]) == np.sign(mr["b_g2"])
mr["dir_consistent"] = mr["dir_consistent"].fillna(False)

# ---------------- 3) BH-FDR ----------------
def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    qv = ranked * n / (np.arange(1, n + 1))
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(qv, 0, 1)
    q[ok] = out
    return q

mr["fdr_fin_all"] = bh(mr["p_fin"].values)   # 全 eGene 参考 (供溯源)

# ---------------- 4) COLOC (近似 coloc.abf / Wakefield ABF) ----------------
def _logsumexp(a):
    a = np.asarray(a, float)
    mx = np.max(a)
    return mx + math.log(np.sum(np.exp(a - mx)))

def _logdiffexp(a, b):
    # log(exp(a) - exp(b)), a>b
    return a + math.log1p(-math.exp(b - a))

def coloc_abf_single(b_e, se_e, b_g, se_g, maf_e, maf_g, N_e, N_g,
                     N_case, p1=1e-4, p2=1e-4, p12=1e-5):
    """
    近似 coloc.abf (independent-SNP approximation):
      ABF_i^e / ABF_i^g = Wakefield ABF.
      H0=1; H1=p1*ΣABF_e; H2=p2*ΣABF_g;
      H3=p1*p2*[(ΣABF_e)(ΣABF_g) - Σ(ABF_e*ABF_g)];
      H4=p12*Σ(ABF_e*ABF_g).
    返回 PPH4, PPH3.
    """
    W = 0.15 ** 2
    def log_abf_single(b, V):
        if V <= 0 or np.isnan(b) or np.isnan(V):
            return 0.0
        z = b / math.sqrt(V)
        r = W / (W + V)
        return 0.5 * math.log(1 - r) + 0.5 * r * z * z
    lbf_e = np.array([log_abf_single(b, s ** 2) for b, s in zip(b_e, se_e)])
    lbf_g = np.array([log_abf_single(b, s ** 2) for b, s in zip(b_g, se_g)])
    le = _logsumexp(lbf_e)                 # log Σ ABF_e
    lg = _logsumexp(lbf_g)                 # log Σ ABF_g
    l_eg = _logsumexp(lbf_e + lbf_g)       # log Σ ABF_e*ABF_g  (H4 core)
    # H3 core = log[(ΣABF_e)(ΣABF_g) - ΣABF_e*ABF_g]
    l3a = le + lg
    l3 = _logdiffexp(l3a, l_eg) if l3a > l_eg else -1e30
    lp1, lp2, lp12 = math.log(p1), math.log(p2), math.log(p12)
    log_post = np.array([0.0, lp1 + le, lp2 + lg, lp1 + lp2 + l3, lp12 + l_eg])
    log_post -= log_post.max()
    post = np.exp(log_post)
    post /= post.sum()
    return post[4], post[3]  # PPH4, PPH3

# 加载可信集 (每基因取其实例 cs_id 的全部变异 beta/se)
cred = pd.read_csv(CRED, sep="\t")
cred["ensg"] = cred["molecular_trait_id"]
cred["rsid"] = cred["rsid"].astype(str)

# 需要 FinnGen 区域 β 的变异集合 = 所有候选基因可信集变异 rsid
need_rsids = set(cred["rsid"])
print(f"[coloc] credible-set variants needed: {len(need_rsids)}", flush=True)

# 扫描 FinnGen 全量, 收集这些 rsid 的 β/se/af
fin_reg = {}
with gzip.open(FIN_FULL, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[4] in need_rsids:
            fin_reg[p[4]] = (float(p[8]), float(p[9]), float(p[6]), p[2], p[3])  # b, se, af_alt, ref, alt
print(f"[coloc] FinnGen regional variants matched: {len(fin_reg)}", flush=True)

PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
coloc_rows = []
for ensg, g in cred.groupby("ensg"):
    sub = g.copy()
    b_e, se_e, b_g, se_g = [], [], [], []
    for _, r in sub.iterrows():
        fr = fin_reg.get(r["rsid"])
        if fr is None:
            continue
        bF, seF, afF, ref, alt = fr
        # FinnGen effect allele = alt; 可信集 variant: chr_pos_ref_alt -> ea=alt
        ea_e = str(r["variant"]).split("_")[3]
        oa_e = str(r["variant"]).split("_")[2]
        if alt == ea_e and ref == oa_e:
            sgn = 1.0
        elif alt == oa_e and ref == ea_e:
            sgn = -1.0
        else:
            continue
        if (ea_e, oa_e) in PAL and 0.4 < afF < 0.6:
            continue
        b_e.append(float(r["beta"])); se_e.append(float(r["se"]))
        b_g.append(bF * sgn); se_g.append(seF)
    if len(b_e) == 0:
        continue
    pph4, pph3 = coloc_abf_single(np.array(b_e), np.array(se_e),
                                  np.array(b_g), np.array(se_g),
                                  None, None, N_EQTL, N_FIN, N_FIN_CASE)
    coloc_rows.append(dict(ensg=ensg, n_snp_coloc=len(b_e), PPH4=pph4, PPH3=pph3))
coloc = pd.DataFrame(coloc_rows)
print(f"[coloc] genes with coloc: {len(coloc)}", flush=True)

mr = mr.merge(coloc, on="ensg", how="left")

def tier(pph4):
    if pd.isna(pph4):
        return "NA"
    if pph4 > 0.80:
        return "Tier1"
    if pph4 > 0.50:
        return "Tier2"
    return "Excluded"
mr["Tier"] = mr["PPH4"].apply(tier)

# ---------------- 4b) 风险优先候选宇宙 (Layer2 溯源: Figure3 MAGMA) ----------------
# MR 候选基因应限于具 PE 基因级遗传关联证据的胎盘 cis-eGene;
# 以 Figure3 fig3A 的 MAGMA 基因级 z (z>=1.0 作为优先筛选, 非显著性判定) 界定 discovery 宇宙。
PRIO = "/workspace/results/figure3/fig3A_gene_priority.csv"
prio = pd.read_csv(PRIO)[["gene", "magma_z", "scdrs_corr", "priority"]]
mr = mr.merge(prio, left_on="symbol", right_on="gene", how="left").drop(columns=["gene"])
mr["is_candidate"] = mr["magma_z"] >= 1.0
mr.loc[mr["magma_z"].isna(), "is_candidate"] = False
mr["fdr_fin"] = np.nan
_m = mr["is_candidate"].values
mr.loc[mr["is_candidate"], "fdr_fin"] = bh(mr.loc[mr["is_candidate"], "p_fin"].values)

mr.to_csv(f"{OUT}/fig4_mr_gene.csv", index=False)
cand_df = mr[mr["is_candidate"]].copy()
cand_df.to_csv(f"{OUT}/fig4_candidate_genes.csv", index=False)
print(f"[universe] 风险优先候选 cis-eGene (MAGMA z>=1.0): {len(cand_df)} / {len(mr)}", flush=True)

# ---------------- 5) 三级漏斗 ----------------
u = cand_df
L0 = len(u)
mac = (u["fdr_fin"] < 0.05)
L1 = int(mac.sum())
mcac = mac & u["dir_consistent"] & (u["p_meta"] < 0.05)
L2 = int(mcac.sum())
L3_t1 = int((mcac & (u["Tier"] == "Tier1")).sum())
L3_t2 = int((mcac & (u["Tier"] == "Tier2")).sum())
excl = int((mcac & ~(u["Tier"].isin(["Tier1", "Tier2"]))).sum())

funnel = pd.DataFrame([
    dict(level="L0: 风险优先候选 cis-eGene", n=L0,
         detail="胎盘 eQTL 可信集 IV (F>10,p<1e-5) ∩ MAGMA z>=1"),
    dict(level="L1: Discovery-MR (FDR<0.05)", n=L1, detail="胎盘 eQTL -> FinnGen R12, BH-FDR<0.05"),
    dict(level="L2: dual-GWAS meta robust", n=L2, detail="方向一致 & 随机效应 meta-p<0.05"),
    dict(level="L3: Tier1 (PPH4>0.80)", n=L3_t1, detail="COLOC 近似 ABF 高置信"),
    dict(level="L3: Tier2 (0.50<PPH4<=0.80)", n=L3_t2, detail="COLOC 中等置信"),
    dict(level="L3: excluded (PPH4<=0.50)", n=excl, detail="共定位不成立, 剔除"),
])
funnel.to_csv(f"{OUT}/fig4_funnel.csv", index=False)

print("\n===== FUNNEL (风险优先候选宇宙) =====")
print(funnel.to_string(index=False))
print("\n===== L1 PASS (FDR<0.05) =====")
print(u[mac].sort_values("p_fin")[["symbol", "b_fin", "p_fin", "fdr_fin", "b_g2", "p_meta",
                                   "dir_consistent", "PPH4", "Tier"]].to_string(index=False))
print(f"\n[out] fig4_mr_gene.csv ({len(mr)}), fig4_candidate_genes.csv ({len(cand_df)}), "
      f"fig4_funnel.csv", flush=True)
