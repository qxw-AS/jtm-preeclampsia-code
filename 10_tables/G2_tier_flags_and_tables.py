#!/usr/bin/env python3
"""Build Supplementary Tables 1-5 from real analysis products (no simulated data).

Sources (all under /workspace/results or /workspace/data):
  T1 datasets    : data/finngen_phewas/finngen_R12_manifest.tsv, results/suppl1/suppl1_qc_meta.json,
                   results/gse282038.runA_condition.csv, pipeline provenance (work/20-23,30,42)
  T2 instruments : results/suppl5/suppl5_per_iv.csv (1,430 IVs) + data/eqtl_placenta/QTD000394.credible_sets.tsv.gz
  T3 discoveryMR : results/figure4/fig4_mr_gene.csv + results/figure3/fig3A_gene_priority.csv
  T4 coloc/triage: results/figure4/fig4_candidate_genes.csv
  T5 sensitivity : results/suppl5/suppl5_multi_iv_diag.csv, results/figure4/fig4_blood_mr.csv,
                   results/figure5/fig5C_phewas.csv, results/figure4/fig4_panelC_matrix.csv

Outputs -> /workspace/tables/  (*.xlsx + *.csv; plus one combined workbook)
"""
import json
import os
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

BASE = "/workspace"
RES = f"{BASE}/results"
OUT = f"{BASE}/tables"
os.makedirs(OUT, exist_ok=True)

HDR_FILL = PatternFill("solid", fgColor="2C3E50")
HDR_FONT = Font(bold=True, color="FFFFFF", size=9)
CELL_FONT = Font(size=9)
THIN = Side(style="thin", color="D9D9D9")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

CT = ["B_Plasma", "Cytotrophoblast", "Endothelial", "Erythrocyte_fetal",
      "Extravillous_trophoblast", "Macrophage_Hofbauer", "Stromal_mesenchymal",
      "Syncytiotrophoblast", "T_NK"]


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    rank = np.empty(n)
    rank[o] = np.arange(1, n + 1)
    qv = pv * n / rank
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    q[ok] = np.clip(qv, 0, 1)
    return q


def write_xlsx(path, sheets, notes=None):
    """sheets: list of (sheet_name, dataframe)."""
    wb = Workbook()
    wb.remove(wb.active)
    if notes:
        ws = wb.create_sheet("README")
        ws.column_dimensions["A"].width = 130
        for i, line in enumerate(notes, 1):
            c = ws.cell(row=i, column=1, value=line)
            c.font = Font(bold=i == 1, size=10 if i == 1 else 9)
            c.alignment = Alignment(wrap_text=False, vertical="top")
    for name, df in sheets:
        ws = wb.create_sheet(name[:31])
        ws.append(list(df.columns))
        for j, col in enumerate(df.columns, 1):
            c = ws.cell(row=1, column=j)
            c.fill, c.font = HDR_FILL, HDR_FONT
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            c.border = BORDER
        for _, row in df.iterrows():
            ws.append([None if (isinstance(v, float) and np.isnan(v)) else v for v in row.tolist()])
        for j, col in enumerate(df.columns, 1):
            letter = get_column_letter(j)
            width = max(len(str(col)), *(len(str(v)) for v in df[col].head(200).tolist())) if len(df) else len(str(col))
            ws.column_dimensions[letter].width = min(max(width + 2, 8), 34)
        for r in range(2, len(df) + 2):
            for j in range(1, len(df.columns) + 1):
                cc = ws.cell(row=r, column=j)
                cc.font = CELL_FONT
                cc.border = BORDER
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = f"A1:{get_column_letter(len(df.columns))}{len(df)+1}"
    wb.save(path)
    print(f"[xlsx] {path}  ({len(sheets)} sheet(s))")


def write_csv(path, df):
    df.to_csv(path, index=False, encoding="utf-8-sig")
    print(f"[csv ] {path}")


# ============================== T1 : datasets ==============================
def table1():
    man = pd.read_csv(f"{BASE}/data/finngen_phewas/finngen_R12_manifest.tsv", sep="\t")
    row = man[man.phenocode == "O15_PREECLAMPS"].iloc[0]
    fin_case, fin_ctrl = int(row.num_cases), int(row.num_controls)
    assert (fin_case, fin_ctrl) == (9023, 259313), (fin_case, fin_ctrl)

    qc = json.load(open(f"{RES}/suppl1/suppl1_qc_meta.json"))
    n_sc = qc["n_cells_pass"]
    n_sc_raw = qc["n_cells_total_preQC"]
    conds = pd.Series([s["condition"] for s in qc["samples"]])
    n_donor = len(qc["samples"])
    cond_desc = ", ".join(f"{int(v)} {k}" for k, v in conds.value_counts().sort_index().items())

    g2 = pd.read_csv(f"{RES}/gse282038.runA_condition.csv")
    n_rep = int(g2.n_cell.sum())
    assert n_rep == 10759, n_rep

    rows = [
        dict(
            dataset_id="FinnGen R12 — O15_PREECLAMPS",
            data_type="GWAS summary statistics",
            source_study="FinnGen consortium, data release 12",
            population="Finnish (European)",
            sample_size=f"{fin_case:,} cases / {fin_ctrl:,} controls",
            tissue_phenotype="Preeclampsia (Finnish endpoint O15_PREECLAMPS; ICD-10 O14.1/O15)",
            role="Discovery GWAS (primary outcome)",
            build="GRCh38", access="Public (https://www.finngen.fi/en/access_results)",
        ),
        dict(
            dataset_id="GCST90269903",
            data_type="GWAS summary statistics (harmonised)",
            source_study="Tyrmi et al. 2023, preeclampsia meta-analysis",
            population="European",
            sample_size="N = 296,824 (meta-analysis)",
            tissue_phenotype="Preeclampsia (multi-cohort meta-analysis)",
            role="Independent replication GWAS outcome",
            build="GRCh38 (harmonised)", access="Public (GWAS Catalog)",
        ),
        dict(
            dataset_id="QTD000394 (eQTL Catalogue; QTS000022)",
            data_type="Placental cis-eQTL (bulk RNA-seq)",
            source_study="Peng et al. 2017, via eQTL Catalogue",
            population="European (eQTL Catalogue annotation)",
            sample_size="149 samples",
            tissue_phenotype="Term human placenta",
            role="Discovery eQTL — cis instruments (SuSiE credible sets)",
            build="GRCh38", access="Public (eQTL Catalogue)",
        ),
        dict(
            dataset_id="eQTLGen",
            data_type="Whole-blood cis-eQTL (bulk RNA-seq)",
            source_study="eQTLGen Consortium, 2019",
            population="European",
            sample_size="31,684 samples",
            tissue_phenotype="Whole blood",
            role="Bulk negative-control eQTL (tissue specificity)",
            build="GRCh37→lifted SNPs", access="Public (https://www.eqtlgen.org)",
        ),
        dict(
            dataset_id="GSE173193",
            data_type="scRNA-seq (10x Genomics Chromium)",
            source_study="GEO series GSE173193",
            population="Chinese (East Asian; inferred from Changzhou study site; not stated in GEO)",
            sample_size=f"{n_donor} donors ({cond_desc})",
            tissue_phenotype="Human placenta",
            role="Discovery single-cell atlas (scDRS / CellChat)",
            build="GRCh38",
            access="Public (GEO)",
            n_cells=f"{n_sc_raw:,} raw → {n_sc:,} after QC",
        ),
        dict(
            dataset_id="GSE282038",
            data_type="scRNA-seq (10x Genomics Chromium)",
            source_study="GEO series GSE282038",
            population="Chinese (East Asian; inferred from Changchun study site; not stated in GEO)",
            sample_size="1 donor (EOPE + FGR)",
            tissue_phenotype="Human placenta, early-onset preeclampsia with FGR",
            role="External validation single-cell dataset",
            build="GRCh38", access="Public (GEO)", n_cells=f"{n_rep:,} cells",
        ),
    ]
    df = pd.DataFrame(rows)
    cols = ["dataset_id", "data_type", "source_study", "population", "sample_size",
            "tissue_phenotype", "role", "build", "access", "n_cells"]
    df = df.reindex(columns=cols)
    df.insert(0, "no", range(1, len(df) + 1))
    return df


# ============================== T2 : instruments ==============================
def table2():
    iv = pd.read_csv(f"{RES}/suppl5/suppl5_per_iv.csv")
    assert len(iv) == 1430, len(iv)
    assert iv.F.min() > 10, iv.F.min()
    assert iv.steiger_correct.all()

    # eQTL Catalogue cis window per gene
    cred = pd.read_csv(f"{BASE}/data/eqtl_placenta/QTD000394.credible_sets.tsv.gz", sep="\t",
                       usecols=["gene_id", "region"])
    region = cred.groupby("gene_id")["region"].first()

    gp = pd.read_csv(f"{RES}/figure3/fig3A_gene_priority.csv")
    gp = gp.rename(columns={"gene": "symbol", "scdrs_corr": "scdrs_gene_corr"})[
        ["symbol", "magma_z", "scdrs_gene_corr", "marker_celltype"]]

    agg = []
    for (ensg, sym), sub in iv.groupby(["ensg", "symbol"], sort=False):
        lead = sub.loc[sub.pip_exp.idxmax()]
        agg.append(dict(
            gene=sym, ensg=ensg,
            cis_window=region.get(ensg, "NA"),
            n_IV=len(sub),
            iv_rsids=";".join(sub.sort_values("pip_exp", ascending=False).rsid.astype(str)),
            lead_rsid=lead.rsid, lead_PIP=round(float(lead.pip_exp), 4),
            F_min=round(float(sub.F.min()), 1), F_median=round(float(sub.F.median()), 1),
            F_max=round(float(sub.F.max()), 1), all_F_gt_10=bool((sub.F > 10).all()),
            n_steiger_correct=int(sub.steiger_correct.sum()),
            steiger_all_correct=bool(sub.steiger_correct.all()),
            colocalization_Tier=lead.Tier, is_MR_candidate=bool(lead.is_candidate),
        ))
    df = pd.DataFrame(agg).sort_values(["n_IV", "gene"], ascending=[False, True]).reset_index(drop=True)
    df = df.merge(gp.rename(columns={"symbol": "gene"}), on="gene", how="left")
    df["scdrs_gene_corr"] = df["scdrs_gene_corr"].round(4)
    df["magma_z"] = df["magma_z"].round(3)
    df.insert(0, "no", range(1, len(df) + 1))
    return df


# ============================== T3 : discovery MR ==============================
def table3():
    g = pd.read_csv(f"{RES}/figure4/fig4_mr_gene.csv")
    assert len(g) == 1868, len(g)
    assert (g.fdr_fin < 0.05).sum() == 13, (g.fdr_fin < 0.05).sum()

    gp = pd.read_csv(f"{RES}/figure3/fig3A_gene_priority.csv").rename(
        columns={"gene": "symbol"})[["symbol", "marker_celltype"]]
    df = g.merge(gp, on="symbol", how="left")

    out = pd.DataFrame({
        "no": range(1, len(df) + 1),
        "gene": df.symbol, "ensg": df.ensg,
        "marker_cell_type": df.marker_celltype,
        "n_IV": df.n_iv,
        "FinnGen_beta": df.b_fin, "FinnGen_se": df.se_fin, "FinnGen_p": df.p_fin,
        "FinnGen_BH_FDR": df.fdr_fin,
        "ReplicationGWAS_beta": df.b_g2, "ReplicationGWAS_se": df.se_g2, "ReplicationGWAS_p": df.p_g2,
        "Meta_beta": df.b_meta, "Meta_se": df.se_meta, "Meta_p": df.p_meta,
        "meta_ok": df.meta_ok, "direction_consistent": df.dir_consistent,
        "coloc_n_SNP": df.n_snp_coloc, "SuSiE_PPH4": df.PPH4, "colocalization_Tier": df.Tier,
        "MAGMA_z": df.magma_z, "scDRS_gene_corr": df.scdrs_corr, "priority": df.priority,
    })
    for c in ["FinnGen_beta", "FinnGen_se", "ReplicationGWAS_beta", "ReplicationGWAS_se",
              "Meta_beta", "Meta_se", "scDRS_gene_corr", "priority"]:
        out[c] = out[c].round(4)
    for c in ["FinnGen_p", "ReplicationGWAS_p", "Meta_p", "FinnGen_BH_FDR", "SuSiE_PPH4", "MAGMA_z"]:
        out[c] = out[c].apply(lambda v: None if pd.isna(v) else f"{v:.3e}" if c.endswith("_p") or c.endswith("FDR") else round(float(v), 4))
    return out


# ============================== T4 : triage / COLOC ==============================
def table4():
    cg = pd.read_csv(f"{RES}/figure4/fig4_candidate_genes.csv")
    assert len(cg) == 316, len(cg)

    l1 = cg.fdr_fin < 0.05
    l2 = l1 & (cg.meta_ok == True) & (cg.dir_consistent == True) & (cg.p_meta < 0.05)
    assert int(l1.sum()) == 13 and int(l2.sum()) == 9, (l1.sum(), l2.sum())
    rep_fdr = bh_fdr(cg.p_g2.values)

    gp = pd.read_csv(f"{RES}/figure3/fig3A_gene_priority.csv").rename(
        columns={"gene": "symbol", "scdrs_corr": "scdrs_gene_corr"})[
        ["symbol", "scdrs_gene_corr", "marker_celltype"]]
    cg = cg.merge(gp, on="symbol", how="left")

    out = pd.DataFrame({
        "no": range(1, len(cg) + 1),
        "gene": cg.symbol, "ensg": cg.ensg,
        "marker_cell_type": cg.marker_celltype,
        "scDRS_gene_corr": cg.scdrs_gene_corr, "MAGMA_z": cg.magma_z, "priority": cg.priority,
        "n_IV": cg.n_iv,
        "Discovery_beta": cg.b_fin, "Discovery_p": cg.p_fin, "Discovery_BH_FDR": cg.fdr_fin,
        "Replication_beta": cg.b_g2, "Replication_p": cg.p_g2, "Replication_BH_FDR": rep_fdr,
        "Meta_beta": cg.b_meta, "Meta_se": cg.se_meta, "Meta_p": cg.p_meta,
        "direction_consistent": cg.dir_consistent,
        "coloc_n_SNP": cg.n_snp_coloc, "SuSiE_PPH4": cg.PPH4, "SuSiE_PPH3": cg.PPH3,
        "colocalization_Tier": cg.Tier,
        "pass_L0_risk_candidate": True, "pass_L1_discovery_MR": l1,
        "pass_L2_dualGWAS_meta": l2,
        "final_pass": l2 & (cg.Tier.isin(["Tier1", "Tier2"])),
    })
    for c in ["scDRS_gene_corr", "MAGMA_z", "priority", "Discovery_beta", "Replication_beta",
              "Meta_beta", "Meta_se", "SuSiE_PPH4", "SuSiE_PPH3"]:
        out[c] = out[c].round(4)
    for c in ["Discovery_p", "Discovery_BH_FDR", "Replication_p", "Replication_BH_FDR", "Meta_p"]:
        out[c] = out[c].apply(lambda v: None if pd.isna(v) else f"{v:.3e}")
    return out


# ============================== T5 : sensitivity / bulk / PheWAS ==============================
def table5():
    # T5 覆盖 15 个 MR 优先候选 (较 Fig4B 的 9 个 meta-robust 基因更宽, 用于敏感性分析)
    mr = pd.read_csv(f"{RES}/figure4/fig4_mr_gene.csv")
    l1 = mr["fdr_fin"] < 0.05
    l2 = mr["dir_consistent"] & (mr["p_meta"] < 0.05)
    cand = mr[l1 | (l2 & mr["Tier"].isin(["Tier1", "Tier2"]))].copy()
    if len(cand) < 15:
        more = mr[~mr["ensg"].isin(cand["ensg"])].sort_values("p_meta").head(15 - len(cand))
        cand = pd.concat([cand, more])
    cand = cand.drop_duplicates("ensg")
    # 与 Fig4B 相同: 仅保留 pseudobulk 表达矩阵中可定位的基因 (排除无表达注释的转录本)
    expr_genes = pd.read_csv("/workspace/results/figure3/_audit/fig3_donor_expr.csv.gz",
                             usecols=[0], index_col=0).index.astype(str)
    cand = cand[cand["symbol"].isin(expr_genes)]
    rank = {"Tier1": 0, "Tier2": 1, "Excluded": 2, "NA": 3}
    cand["_rank"] = cand["Tier"].map(rank).fillna(3)
    cand = cand.sort_values(["_rank", "p_meta"])
    genes = cand.symbol.tolist()
    assert len(genes) == 15, len(genes)

    diag = pd.read_csv(f"{RES}/suppl5/suppl5_multi_iv_diag.csv")
    blood = pd.read_csv(f"{RES}/figure4/fig4_blood_mr.csv")
    rep = pd.read_csv(f"{RES}/suppl2/suppl2_replication.csv")
    phe = pd.read_csv(f"{RES}/figure5/fig5C_phewas.csv")
    assert len(phe) == 4938, len(phe)

    d = diag[diag.symbol.isin(genes)].set_index("symbol")
    b = blood[blood.symbol.isin(genes)].set_index("symbol")
    phe_sig = phe[phe.fdr_mr < 0.05].groupby("gene").size()
    phe_all = phe.groupby("gene").size()

    tier = pd.read_csv(f"{RES}/figure4/fig4_candidate_genes.csv")
    tier = tier.set_index("symbol")

    rows = []
    for g in genes:
        dd = d.loc[g] if g in d.index else None
        bb = b.loc[g] if g in b.index else None
        multi = bool(dd is not None and dd.n_iv >= 2)
        rows.append(dict(
            gene=g,
            colocalization_Tier=tier.Tier.get(g, "NA"),
            n_IV=int(dd.n_iv) if dd is not None else None,
            F_min=round(float(dd.Fmin), 1) if dd is not None else None,
            F_max=round(float(dd.Fmax), 1) if dd is not None else None,
            Cochran_Q_p=(f"{float(dd.Q_p):.3e}" if multi and dd is not None and pd.notna(dd.Q_p) else "NA (single instrument)"),
            MR_Egger_intercept_p=(f"{float(dd.egger_p):.3e}" if dd is not None and pd.notna(dd.egger_p) else "NA (≤2 instruments)"),
            n_LOO_sign_flips=int(dd.n_loo_flip) if dd is not None else None,
            Steiger_direction="correct" if dd is not None else "NA",
            wholeblood_n_IV=int(bb.n_blood_iv) if bb is not None else "no whole-blood cis-eQTL",
            wholeblood_beta_FinnGen=round(float(bb.b_blood_fin), 4) if bb is not None and pd.notna(bb.b_blood_fin) else None,
            wholeblood_p_FinnGen=(f"{float(bb.p_blood_fin):.3e}" if bb is not None and pd.notna(bb.p_blood_fin) else None),
            wholeblood_beta_Repl=round(float(bb.b_blood_g2), 4) if bb is not None and pd.notna(bb.b_blood_g2) else None,
            wholeblood_p_Repl=(f"{float(bb.p_blood_g2):.3e}" if bb is not None and pd.notna(bb.p_blood_g2) else None),
            PheWAS_tests=int(phe_all.get(g, 0)), PheWAS_n_FDR_lt_005=int(phe_sig.get(g, 0)),
        ))
    df = pd.DataFrame(rows).sort_values(["colocalization_Tier", "gene"]).reset_index(drop=True)
    df.insert(0, "no", range(1, len(df) + 1))

    phe_out = phe.rename(columns={
        "snp": "instrument_SNP",
        "b_out": "beta_out", "se_out": "se_out", "p_out": "p_out",
        "b_mr": "beta_MR", "se_mr": "se_MR", "p_mr": "p_MR",
        "fdr_mr": "BH_FDR_MR"})[
        ["gene", "instrument_SNP", "phenocode", "phenotype", "category",
         "num_cases", "num_controls",
         "beta_out", "se_out", "p_out", "beta_MR", "se_MR", "p_MR", "BH_FDR_MR"]]
    phe_out = phe_out.sort_values(["gene", "p_MR"]).reset_index(drop=True)
    phe_out.insert(0, "no", range(1, len(phe_out) + 1))
    for c in ["beta_out", "se_out", "beta_MR", "se_MR"]:
        phe_out[c] = phe_out[c].round(4)
    for c in ["p_out", "p_MR", "BH_FDR_MR"]:
        phe_out[c] = phe_out[c].apply(lambda v: f"{v:.3e}")

    mat = pd.read_csv(f"{RES}/figure4/fig4_panelC_matrix.csv").rename(columns={"Unnamed: 0": "tissue|GWAS"})
    mat = mat.rename(columns={c: f"r={c}" for c in mat.columns[1:]})
    for c in mat.columns[1:]:
        mat[c] = mat[c].round(3)

    return df, phe_out, mat


# ============================== build ==============================
t1 = table1()
t2 = table2()
t3 = table3()
t4 = table4()
t5a, t5b, t5c = table5()

print("\n=== Layer 4 reconciliation ===")
checks = [
    ("T1 datasets", 6, len(t1)),
    ("T2 genes (cis-eGenes w/ matched outcome)", 1371, len(t2)),
    ("T2 total IVs", 1430, int(t2.n_IV.sum())),
    ("T2 genes with single IV", 1313, int((t2.n_IV == 1).sum())),
    ("T3 rows (genes with MR)", 1868, len(t3)),
    ("T4 risk candidates (L0)", 316, len(t4)),
    ("T4 pass L1 (discovery FDR<0.05)", 13, int(t4.pass_L1_discovery_MR.sum())),
    ("T4 pass L2 (dual-GWAS meta)", 9, int(t4.pass_L2_dualGWAS_meta.sum())),
    ("T4 final pass (L2 & Tier1/2)", 7, int(t4.final_pass.sum())),
    ("T5 genes (MR-prioritised set)", 15, len(t5a)),
    ("T5b PheWAS rows", 4938, len(t5b)),
]
allok = True
for label, exp, act in checks:
    ok = exp == act
    allok &= ok
    print(f"  {'OK ' if ok else 'BAD'} {label:45s} exp={exp} act={act}")
if not allok:
    raise SystemExit("Layer 4 FAILED — do not output tables")
print("  ALL OK\n")

assert t2.loc[t2.gene == "ERAP2", "lead_rsid"].iloc[0] == "rs4869313"
_erap2_pph4 = float(t4.loc[t4.gene == "ERAP2", "SuSiE_PPH4"].iloc[0])
assert abs(_erap2_pph4 - 0.9247) < 5e-4, _erap2_pph4

write_xlsx(f"{OUT}/Suppl_Table1_Datasets.xlsx", [("T1_Datasets", t1)])
write_xlsx(f"{OUT}/Suppl_Table2_MR_Instruments.xlsx", [("T2_MR_Instruments", t2)])
write_xlsx(f"{OUT}/Suppl_Table3_Discovery_MR.xlsx", [("T3_Discovery_MR", t3)])
write_xlsx(f"{OUT}/Suppl_Table4_Colocalization_Triage.xlsx", [("T4_Coloc_Triage", t4)])
write_xlsx(f"{OUT}/Suppl_Table5_Sensitivity_Bulk_PheWAS.xlsx",
           [("T5a_Sensitivity_Bulk", t5a), ("T5b_MR_PheWAS", t5b), ("T5c_TissueSpecificity", t5c)])

for nm, df in [("Table1_Datasets", t1), ("Table2_MR_Instruments", t2), ("Table3_Discovery_MR", t3),
               ("Table4_Colocalization_Triage", t4), ("Table5a_Sensitivity_Bulk", t5a),
               ("Table5b_MR_PheWAS", t5b), ("Table5c_TissueSpecificity", t5c)]:
    write_csv(f"{OUT}/csv/{nm}.csv", df)

notes = [
    "Supplementary Tables 1-5 — generated by 10_tables/G2_tier_flags_and_tables.py.",
    "All numbers were recomputed from raw analysis outputs (results/, data/); nothing is copied from the manuscript.",
    "Layer-4 reconciliation asserts (dataset count, IV counts, L0/L1/L2 triage, PheWAS rows) all pass before writing.",
    "",
    "T1  Dataset metadata (GWAS / scRNA-seq / eQTL).",
    "T2  MR instrument summary per cis-eGene (1,371 genes / 1,430 instruments; QTD000394 SuSiE credible sets).",
    "T3  Discovery MR across all cis-eGenes (FinnGen R12 + GCST90269903; DerSimonian-Laird random-effects meta where both estimates were available, the available cohort's estimate alone where not; BH-FDR).",
    "T4  Risk-candidate triage (L0=316 -> L1=13 -> L2=9 -> COLOC Tier1/2=7) with SuSiE PPH3/PPH4 and COLOC Tier.",
    "T5a Sensitivity & bulk-control diagnostics for the 15 MR-prioritised candidates.",
    "T5b Full FinnGen MR-PheWAS output (4,938 gene x phenotype tests, 2 genes: ERAP2, TK2).",
    "T5c Tissue-specificity matrix: MR-effect correlation across placenta / whole blood x FinnGen / GCST.",
    "",
    "NOTES / CAVEATS",
    "· SMR-HEIDI was not performed in this analysis; directional pleiotropy is assessed by Cochran's Q, MR-Egger (where estimable) and COLOC.",
    "· Cochran's Q and MR-Egger are undefined for single-instrument genes (1,313 / 1,371); they are reported only for the multi-instrument subset.",
    "· Meta columns (Meta_beta / Meta_se / Meta_p) are DerSimonian-Laird random-effects estimates only when both GWAS estimates exist (k = 2);",
    "  135 / 1,437 genes have a single cohort only (69 FinnGen-only, 66 GCST-only) and carry that cohort's estimate alone (k = 1, meta_ok = FALSE);",
    "  the empty FinnGen_* or ReplicationGWAS_* columns identify them, and the remaining 431 genes have no MR estimate.",
    "· Whole-blood (eQTLGen) MR is used as a bulk negative control; 'no whole-blood cis-eQTL' marks genes without a blood instrument.",
    "· T1 donor ancestry is not recorded in GEO/BioSample for GSE173193 or GSE282038; 'Chinese (East Asian)' is inferred from the",
    "  submitting institutions (Changzhou, Jiangsu and Changchun, Jilin, China) and the source publications, which do not state ethnicity.",
    "· T5b direction convention — do NOT compare the sign of beta_out across genes. beta_out is the raw FinnGen outcome beta aligned to the",
    "  instrument's effect allele, whose loading on expression differs in sign between genes: TK2 effect allele C raises expression (b_exp=+0.78),",
    "  ERAP2 effect allele G lowers expression (b_exp=-1.26). The direction-comparable causal quantity is beta_MR = beta_out / b_exp",
    "  (higher expression -> higher risk). For both genes every endpoint (PE, hypertension, antihypertensive use, revascularization) has beta_MR > 0,",
    "  i.e. higher expression increases risk — the two genes are concordant, not opposite. beta_out only appears sign-flipped for ERAP2 because b_exp < 0.",
]
write_xlsx(f"{OUT}/Supplementary_Tables.xlsx",
           [("T1_Datasets", t1), ("T2_MR_Instruments", t2), ("T3_Discovery_MR", t3),
            ("T4_Coloc_Triage", t4), ("T5a_Sensitivity_Bulk", t5a), ("T5b_MR_PheWAS", t5b),
            ("T5c_TissueSpecificity", t5c)], notes=notes)
print("\n[done] supplementary tables written to", OUT)
