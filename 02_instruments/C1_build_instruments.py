#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure4 上游 (1/6): 构建 MR 工具变量表 (Discovery = Peng_2017 胎盘 cis-eQTL)
=====================================================================
数据源 (Layer 2 溯源):
  · eQTL Catalogue QTS000022 / QTD000394 (Peng_2017 胎盘, ge, n=149, GRCh38)
    SuSiE 可信集: QTD000394.credible_sets.tsv.gz
    列: molecular_trait_id, gene_id, cs_id, variant, rsid, cs_size, pip, pvalue, beta, se, z, cs_min_r2, region
  · GENCODE v44 (GRCh38) 提供 ENSG -> symbol / TSS
工具变量规则:
  · 每个基因取每个可信集(cs_id)内 PIP 最高的变异 (SuSiE 已保证 LD 独立, 替代 LD-clumping)
  · p<1e-5; type=SNP; F=(beta/se)^2 > 10 过滤弱工具
输出: /workspace/results/figure4/fig4_instruments.csv
"""
import os, gzip, re
import numpy as np
import pandas as pd

OUT = "/workspace/results/figure4"
os.makedirs(OUT, exist_ok=True)
CRED = "/workspace/data/eqtl_placenta/QTD000394.credible_sets.tsv.gz"
GTF = "/workspace/data/annotation/gencode.v44.annotation.gtf.gz"

# ---------------- 1) ENSG -> symbol / TSS (GRCh38) ----------------
gene_info = {}
with gzip.open(GTF, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        chrom, start, end, strand, attr = parts[0], int(parts[3]), int(parts[4]), parts[6], parts[8]
        m = re.search(r'gene_id "([^"]+)"', attr)
        n = re.search(r'gene_name "([^"]+)"', attr)
        if not m:
            continue
        gid = m.group(1).split(".")[0]
        gname = n.group(1) if n else gid
        tss = start if strand == "+" else end
        # 主染色体
        if chrom in {f"chr{i}" for i in range(1, 23)}:
            gene_info[gid] = (gname, chrom.replace("chr", ""), tss)
print(f"[gencode] genes mapped: {len(gene_info)}")

# ---------------- 2) SuSiE 可信集 -> 工具变量 ----------------
cred = pd.read_csv(CRED, sep="\t")
cred["gene"] = cred["molecular_trait_id"]
n_genes_all = cred["gene"].nunique()
n_cs_all = cred["cs_id"].nunique()

# 每 cs_id 取 PIP 最高变异
cred = cred.sort_values("pip", ascending=False)
inst = cred.groupby("cs_id", as_index=False).first()

# 过滤
inst = inst[inst["pvalue"] < 1e-5].copy()
inst["type"] = inst["variant"].str.count("_").apply(lambda k: "SNP" if k == 3 else "INDEL")
inst = inst[inst["type"] == "SNP"].copy()
inst["F"] = (inst["beta"] / inst["se"]) ** 2
inst = inst[inst["F"] > 10].copy()

# 拆分 variant -> chrom/pos/ref/alt
vc = inst["variant"].str.split("_", expand=True)
inst["chrom"] = vc[0].str.replace("chr", "", regex=False)
inst["pos"] = vc[1].astype(int)
inst["ea"] = vc[3]   # effect allele = alt
inst["oa"] = vc[2]   # other allele = ref

inst["symbol"] = inst["gene"].map(lambda g: gene_info.get(g, (None, None, None))[0])
inst["tss"] = inst["gene"].map(lambda g: gene_info.get(g, (None, None, None))[2])

inst = inst.rename(columns={"molecular_trait_id": "ensg", "cs_id": "cs",
                            "beta": "b_exp", "se": "se_exp", "pvalue": "p_exp",
                            "z": "z_exp", "pip": "pip_exp"})
cols = ["ensg", "symbol", "chrom", "pos", "tss", "variant", "rsid", "cs", "cs_size",
        "pip_exp", "p_exp", "b_exp", "se_exp", "z_exp", "F", "ea", "oa"]
inst = inst[cols].sort_values("ensg").reset_index(drop=True)
inst.to_csv(f"{OUT}/fig4_instruments.csv", index=False)

print(f"[credible] eGenes={n_genes_all}  credible_sets={n_cs_all}")
print(f"[instruments] genes={inst['ensg'].nunique()}  ivs={len(inst)}  "
      f"p<1e-5 & F>10 & SNP")
print(f"[instruments] F range: {inst['F'].min():.1f} - {inst['F'].max():.1f}; "
      f"multi-IV genes={ (inst.groupby('ensg').size()>1).sum() }")
print(f"[instruments] symbol missing={inst['symbol'].isna().sum()}")
print("\nTop by |b_exp|:")
print(inst.reindex(inst["b_exp"].abs().sort_values(ascending=False).index)
      .head(8)[["ensg", "symbol", "rsid", "b_exp", "se_exp", "p_exp", "F"]].to_string(index=False))
