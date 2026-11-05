#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Central path configuration (checklist item A3).

Every absolute path used by the P0 analysis scripts is collected here.
The archived scripts still contain the original absolute paths (/workspace/...)
exactly as they were run; when redeploying, either

  (a) set PE_PROJECT_ROOT and replace the literals in the scripts with the
      constants below, or
  (b) add `from config import ...` at the top of each script.

Directory layout assumed (created by the download / QC steps):

  $PE_PROJECT_ROOT/
    data/
      finngen/            FinnGen R12 endpoint summary stats + MAGMA p-value input
      gwas2/              GCST90269903 harmonised summary stats
      eqtl_placenta/      QTD000394 SuSiE credible sets
      annotation/         GENCODE v44 GTF, NCBI37.3 gene locations
      ref1000g/           MAGMA 1000G EUR reference panel (g1000_eur.{bed,bim,fam})
      magma/              MAGMA binaries' output
      finngen_phewas/     FinnGen R12 endpoint manifest
      gse173193/          discovery placental scRNA-seq (processed)
      gse282038/          external validation sample (processed)
      scdrs/              scDRS score files
    results/              analysis outputs
    work/                 MAGMA binary + intermediate set files
    tables/               supplementary tables
"""
import os

PROJECT_ROOT = os.environ.get("PE_PROJECT_ROOT", "/workspace")

DATA = os.path.join(PROJECT_ROOT, "data")
RESULTS = os.path.join(PROJECT_ROOT, "results")
WORK = os.path.join(PROJECT_ROOT, "work")
TABLES = os.path.join(PROJECT_ROOT, "tables")

# ---------------------------------------------------------------- GWAS / eQTL
FINNGEN_DIR = os.path.join(DATA, "finngen")
FINNGEN_O15 = os.path.join(FINNGEN_DIR, "finngen_R12_O15_PREECLAMPS.gz")
FINNGEN_MAGMA_PVAL = os.path.join(FINNGEN_DIR, "finngen_O15_PREECLAMPS.magma.pval")

GWAS2_DIR = os.path.join(DATA, "gwas2")
GWAS2_HARMONISED = os.path.join(GWAS2_DIR, "GCST90269903.h.tsv.gz")

EQTL_DIR = os.path.join(DATA, "eqtl_placenta")
EQTL_CREDIBLE_SETS = os.path.join(EQTL_DIR, "QTD000394.credible_sets.tsv.gz")

ANNOT_DIR = os.path.join(DATA, "annotation")
GENCODE_GTF = os.path.join(ANNOT_DIR, "gencode.v44.annotation.gtf.gz")
GENE_LOC_HG19 = os.path.join(ANNOT_DIR, "NCBI37.3.gene.loc")

FINNGEN_PHEWAS_DIR = os.path.join(DATA, "finngen_phewas")
FINNGEN_PHEWAS_MANIFEST = os.path.join(FINNGEN_PHEWAS_DIR, "finngen_R12_manifest.tsv")

# ---------------------------------------------------------------- MAGMA
MAGMA_DIR = os.path.join(DATA, "magma")
MAGMA_BIN = os.path.join(WORK, "magma", "magma")
MAGMA_REF = os.path.join(DATA, "ref1000g", "g1000_eur")
MAGMA_GENES_OUT = os.path.join(MAGMA_DIR, "finngen_O15.genes.out")
MAGMA_GENES_RAW = os.path.join(MAGMA_DIR, "finngen_O15.genes.raw")
MAGMA_PVAL_TSV = os.path.join(MAGMA_DIR, "finngen_O15.pval.tsv")
MAGMA_CELLTYPE_OUT = os.path.join(MAGMA_DIR, "finngen_O15.celltype")
MAGMA_GS = os.path.join(MAGMA_DIR, "O15_PREECLAMPS.gs")
MAGMA_CELLTYPE_SETS = os.path.join(WORK, "magma_celltype_sets.txt")

# ---------------------------------------------------------------- single cell
GSE173193_DIR = os.path.join(DATA, "gse173193")
GSE173193_RAW = os.path.join(GSE173193_DIR, "combined_raw.h5ad")
GSE173193_LOGNORM = os.path.join(GSE173193_DIR, "combined_lognorm.h5ad")
GSE173193_HVG = os.path.join(GSE173193_DIR, "combined_hvg.h5ad")
GSE173193_EMB = os.path.join(GSE173193_DIR, "combined_emb.h5ad")
GSE173193_ANNOT = os.path.join(GSE173193_DIR, "combined_annotated.h5ad")
GSE173193_SCDRS = os.path.join(GSE173193_DIR, "combined_scdrs.h5ad")

GSE282038_DIR = os.path.join(DATA, "gse282038")
GSE282038_RAW = os.path.join(GSE282038_DIR, "adata_raw.h5ad")
GSE282038_INTEGRATED = os.path.join(GSE282038_DIR, "adata_integrated.h5ad")
GSE282038_ANNOT = os.path.join(GSE282038_DIR, "adata_annotated.h5ad")

SCDRS_DIR = os.path.join(DATA, "scdrs")
SCDRS_GSE173193 = os.path.join(SCDRS_DIR, "gse173193")
SCDRS_GSE282038 = os.path.join(SCDRS_DIR, "gse282038")
SCDRS_FULL_SCORE = os.path.join(SCDRS_GSE173193, "O15_PREECLAMPS.full_score.gz")

# ---------------------------------------------------------------- results
FIG3_DIR = os.path.join(RESULTS, "figure3")
FIG3_PRIORITY = os.path.join(FIG3_DIR, "fig3A_gene_priority.csv")
FIG3_DONOR_EXPR = os.path.join(FIG3_DIR, "_audit", "fig3_donor_expr.csv.gz")

FIG4_DIR = os.path.join(RESULTS, "figure4")
FIG4_INSTRUMENTS = os.path.join(FIG4_DIR, "fig4_instruments.csv")
FIG4_OUTCOME_FINNGEN = os.path.join(FIG4_DIR, "fig4_outcome_finngen.csv")
FIG4_OUTCOME_GWAS2 = os.path.join(FIG4_DIR, "fig4_outcome_gwas2.csv")
FIG4_MR_GENE = os.path.join(FIG4_DIR, "fig4_mr_gene.csv")
FIG4_CANDIDATES = os.path.join(FIG4_DIR, "fig4_candidate_genes.csv")
FIG4_FUNNEL = os.path.join(FIG4_DIR, "fig4_funnel.csv")
FIG4_BLOOD_MR = os.path.join(FIG4_DIR, "fig4_blood_mr.csv")            # P1 (23)
FIG4_PANELC_MATRIX = os.path.join(FIG4_DIR, "fig4_panelC_matrix.csv")   # P1 (23)

FIG5_DIR = os.path.join(RESULTS, "figure5")
FIG5_PHEWAS = os.path.join(FIG5_DIR, "fig5C_phewas.csv")

SUPPL5_DIR = os.path.join(RESULTS, "suppl5")
SUPPL5_DIAG = os.path.join(SUPPL5_DIR, "suppl5_multi_iv_diag.csv")
SUPPL5_PER_IV = os.path.join(SUPPL5_DIR, "suppl5_per_iv.csv")

RESULTS_GSE173193 = os.path.join(RESULTS, "gse173193")
RESULTS_GSE282038 = os.path.join(RESULTS, "gse282038")
