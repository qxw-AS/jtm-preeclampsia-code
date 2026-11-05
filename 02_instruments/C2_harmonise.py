#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure4 上游 (2/6): 查询工具变量在结局 GWAS 中的效应并 harmonise
  · GWAS-1 FinnGen R12 O15_PREECLAMPS (本地, GRCh38)
  · GWAS-2 GCST90269903 harmonised (本地 tabix, GRCh38)
输出: fig4_outcome_finngen.csv / fig4_outcome_gwas2.csv
"""
import os, gzip, sys
import numpy as np
import pandas as pd
import pysam

OUT = "/workspace/results/figure4"
os.makedirs(OUT, exist_ok=True)
inst = pd.read_csv(f"{OUT}/fig4_instruments.csv")
inst["rsid"] = inst["rsid"].astype(str)
rsids = set(inst["rsid"])
print(f"[in] instruments: {len(inst)} IVs / {inst['ensg'].nunique()} genes", flush=True)

PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

# ---------------- GWAS-1 FinnGen ----------------
rows = []
seen = 0
with gzip.open("/workspace/data/finngen/finngen_R12_O15_PREECLAMPS.gz", "rt") as f:
    f.readline()
    for line in f:
        seen += 1
        if seen % 5000000 == 0:
            print(f"  [FinnGen] scanned {seen/1e6:.0f}M rows, matched {len(rows)}", flush=True)
        p = line.rstrip("\n").split("\t")
        if p[4] not in rsids:
            continue
        rows.append((p[4], p[0], int(p[1]), p[2], p[3], float(p[8]), float(p[9]), float(p[6]), float(p[10])))
fin = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "ref", "alt", "b_out", "se_out", "p_out", "af_alt"])
print(f"[FinnGen] scanned {seen/1e6:.0f}M rows, matched IVs: {len(fin)}", flush=True)

# harmonise (FinnGen: effect allele = alt)
m1 = inst.merge(fin, on="rsid", how="left")
m1 = m1[m1["b_out"].notna()].copy()
b_new, keep_flag = [], []
for _, r in m1.iterrows():
    a, b, c, d = r["ea"], r["oa"], r["alt"], r["ref"]
    if c == a and d == b:
        b_new.append(r["b_out"]); keep_flag.append(True)
    elif c == b and d == a:
        b_new.append(-r["b_out"]); keep_flag.append(True)
    else:
        b_new.append(np.nan); keep_flag.append(False)
m1["b_out_h"] = b_new
m1 = m1[np.array(keep_flag)].copy()
# 回文歧义
pal_amb = m1.apply(lambda r: (r["ea"], r["oa"]) in PAL and 0.4 < r["af_alt"] < 0.6, axis=1)
m1 = m1[~pal_amb].copy()
m1.to_csv(f"{OUT}/fig4_outcome_finngen.csv", index=False)
print(f"[FinnGen] harmonised kept: {len(m1)} (dropped palindromic/ambiguous={pal_amb.sum()})", flush=True)

# ---------------- GWAS-2 GCST90269903 (远程 tabix, 并行 range 请求) ----------------
# 本地整文件 375MB, 代理下载仅 ~29KB/s (~3.7h); 改用 EBI 支持 Range 的 .tbi,
# 仅按 IV 位点并行取窗口, 8~16 线程 -> 每个 IV ~0.4s, 全量约 15 min.
G2_LOCAL = "/workspace/data/gwas2/GCST90269903.h.tsv.gz"
G2_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
          "GCST90269001-GCST90270000/GCST90269903/harmonised/GCST90269903.h.tsv.gz")
print(f"[GCST] local partial={os.path.exists(G2_LOCAL)} -> 使用远程 tabix 并行查询", flush=True)

import threading
from concurrent.futures import ThreadPoolExecutor
_tls = threading.local()

def _tab():
    t = getattr(_tls, "tab", None)
    if t is None:
        t = pysam.TabixFile(G2_URL, index=G2_URL + ".tbi")
        _tls.tab = t
    return t

def _fetch(row):
    ch, pos = str(row["chrom"]), int(row["pos"])
    recs = []
    for attempt in range(3):
        try:
            recs = list(_tab().fetch(ch, pos - 1, pos + 1))
            break
        except Exception:
            _tls.tab = None
    out = []
    for rec in recs:
        c = rec.split("\t")
        if int(c[1]) != pos:
            continue
        # cols: chrom,pos,eff,oth,beta,se,eaf,pval,rsid,...
        out.append((row["ensg"], row["symbol"], row["rsid"], row["cs"], ch, pos,
                    row["ea"], row["oa"], row["b_exp"], row["se_exp"], row["p_exp"], row["F"],
                    c[2], c[3], float(c[4]), float(c[5]), float(c[6]), float(c[7])))
    return out

g2rows = []
with ThreadPoolExecutor(max_workers=16) as ex:
    done = 0
    for res in ex.map(_fetch, [r for _, r in inst.iterrows()]):
        g2rows.extend(res)
        done += 1
        if done % 200 == 0:
            print(f"  [GCST] queried {done}/{len(inst)} IVs, matched records {len(g2rows)}", flush=True)
g2 = pd.DataFrame(g2rows, columns=["ensg", "symbol", "rsid", "cs", "chrom", "pos", "ea", "oa",
                                   "b_exp", "se_exp", "p_exp", "F", "ea_out", "oa_out",
                                   "b_out", "se_out", "eaf_out", "p_out"])
if len(g2):
    g2 = g2.drop_duplicates(subset=["rsid", "ea_out", "oa_out"]).reset_index(drop=True)
print(f"[GCST] matched IV records: {len(g2)}", flush=True)
out = []
for _, r in g2.iterrows():
    a, b = r["ea"], r["oa"]
    if r["ea_out"] == a and r["oa_out"] == b:
        sign = 1.0
    elif r["ea_out"] == b and r["oa_out"] == a:
        sign = -1.0
    else:
        continue
    if (a, b) in PAL and 0.4 < r["eaf_out"] < 0.6:
        continue
    out.append({**r.to_dict(), "b_out": r["b_out"] * sign})
m2 = pd.DataFrame(out).drop(columns=["eaf_out"], errors="ignore")
m2.to_csv(f"{OUT}/fig4_outcome_gwas2.csv", index=False)
print(f"[GCST] harmonised kept: {len(m2)}", flush=True)
print("done.", flush=True)
