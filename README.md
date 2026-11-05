# jtm-preeclampsia-code

Analysis code for:

> **Placental genetic risk mapping associates ERAP2 and TK2 with non-trophoblast cell states and interface remodelling in preeclampsia**
> Xiaowei Qiu, Department of Obstetrics and Gynecology, Nanfang Hospital, Southern Medical University, Guangzhou, China

The repository reproduces the primary numbers of the manuscript: the instrument
and harmonisation cascade (1,981 lead variants / 1,868 *cis*-eGenes → 1,371 genes
/ 1,430 instruments), the discovery MR with FDR control (316 → 13 candidates), the
dual-GWAS random-effects meta-analysis with overlap correction (9 direction-consistent
genes; 64 → 59 at meta-P < 0.05), the colocalisation tiers (6 Tier 1, 9 Tier 2,
7 retained), the single-cell QC / Harmony / annotation of 55,414 placental cells,
scDRS and MAGMA compartment testing, and the FinnGen MR-PheWAS of the two Tier 1 genes.

---

## 1. Layout

Directories are numbered in run order; file prefixes are the checklist identifiers
used in `docs/CODE_MAP.md`, which maps every script to the manuscript section and to
the numbers it produces.

```
.
├── README.md
├── LICENSE                       MIT
├── CITATION.cff
├── environment.yml               conda environment
├── config.py                     all absolute paths in one place
├── docs/
│   └── CODE_MAP.md               script -> manuscript -> expected numbers
├── 01_magma/
│   ├── K1a_magma_input_pval.py       FinnGen R12 -> MAGMA --pval input
│   ├── K1b_magma_gene_z.sh           MAGMA annotate + gene analysis (v1.10)
│   └── K1c_magma_genesout_to_pval.py entrez -> symbol; scDRS p-value file
├── 02_instruments/
│   ├── C1_build_instruments.py       SuSiE credible sets -> instruments
│   └── C2_harmonise.py               outcome lookup + allele harmonisation
├── 03_mr_meta_coloc/
│   └── D1-D2-F1-G1_mr_meta_coloc.py  Wald / IVW, DL meta, coloc ABF, tiers, funnel
├── 04_overlap_correction/
│   └── F2_overlap_correction.py      replication SE x 1.19, tau^2 fixed
├── 05_sensitivity/
│   └── E1_mr_sensitivity.py          Cochran's Q, MR-Egger, leave-one-out, Steiger
├── 06_singlecell/
│   ├── I1a_sc_qc_build_adata.py      per-sample QC, on-disk concatenation
│   ├── I1b_sc_integrate_harmony.py   HVG -> PCA -> Harmony -> Leiden -> UMAP
│   └── I1c_sc_annotate.py            marker scoring, nine compartments
├── 07_scdrs/
│   ├── J1a_scdrs_input.py            raw counts + metadata for scDRS
│   ├── J1b_scdrs_score.py            scDRS cell-level scoring
│   └── J1c_scdrs_downstream.py       compartment MC association, donor aggregation
├── 08_compartments/
│   └── J2_compartment_magma_sets.py  nine marker sets for MAGMA set-annot
├── 09_phewas/
│   └── M1_mr_phewas.py               FinnGen R12 Results API, 2,469 endpoints
└── 10_tables/
    └── G2_tier_flags_and_tables.py   stage flags + supplementary tables 1-5
```

---

## 2. Data

Raw data are not redistributed; all resources are public.

| Resource | Accession / source | Role |
|---|---|---|
| Placental *cis*-eQTL, SuSiE credible sets | eQTL Catalogue **QTD000394** (QTS000022; Peng 2017, n = 149) | exposure instruments |
| Preeclampsia GWAS, discovery | **FinnGen R12** endpoint `O15_PREECLAMPS` (9,023 / 259,313) | discovery outcome |
| Preeclampsia GWAS, replication | **GCST90269903** (Tyrmi 2023, N = 296,824) | replication outcome |
| FinnGen endpoint manifest | FinnGen R12 Results API | 2,469 endpoints for MR-PheWAS |
| Placental scRNA-seq, discovery | **GSE173193** (58,577 cells → 55,414; 8 donors) | scDRS / CellChat |
| Placental scRNA-seq, external check | **GSE282038**, GSM8634701 (1 donor, EOPE + FGR) | one-donor comparison |
| LD reference for MAGMA | 1000G EUR (`g1000_eur`, hg19) | gene analysis / set test |
| Gene models | GENCODE v44 (GRCh38); NCBI37.3 gene locations (hg19) | annotation |

Download with `GEOparse`, `ffq` or the SRA Toolkit for the GEO series; FinnGen R12
summary statistics and the endpoint manifest from the FinnGen release page; GWAS
Catalog harmonised files for GCST90269903.

---

## 3. Installation

```bash
conda env create -f environment.yml
conda activate jtm-preeclampsia
```

Two dependencies are not installable from `environment.yml`:

* **MAGMA v1.10** — official binary, with the 1000G EUR reference panel.
* **CellChat v1.6** — R package; needed only for the cell–cell communication analysis,
  which is not part of this repository (see §6).

The single-cell stack is specified with minimum versions. To freeze exact builds:

```bash
conda env export --no-builds > environment.lock.yml
```

---

## 4. Configuration

Every path used by the analysis is collected in `config.py`. The scripts were run in a
container where the project root was `/workspace`, and the literal paths are preserved
exactly as executed so that the code matches the analysis that produced the manuscript.

To redeploy elsewhere, either mount your data at `/workspace` or rewrite the literals:

```bash
# from the repository root
sed -i 's#/workspace#/your/project/root#g' */*.py */*.sh
```

`config.py` documents the expected directory layout and exposes the same paths as
names, so a script can alternatively `from config import FINNGEN_O15, EQTL_CREDIBLE_SETS`
once the repository root is on `PYTHONPATH`.

---

## 5. Run order

### Step 1 — MAGMA gene-level z-scores and cell-type sets (K1, J2)

```bash
cd 01_magma
python3 K1a_magma_input_pval.py        # FinnGen R12 -> MAGMA --pval input
bash    K1b_magma_gene_z.sh            # annotate window=10,10; gene analysis
python3 K1c_magma_genesout_to_pval.py  # entrez -> symbol; scDRS p-value file

cd ../08_compartments
python3 J2_compartment_magma_sets.py   # nine cell-type marker sets
# competitive gene-set test (see README of that step):
#   magma --gene-results data/magma/finngen_O15.genes.raw \
#         --set-annot work/magma_celltype_sets.txt \
#         --out data/magma/finngen_O15.celltype
```

Edit `MAGMA=` in `K1b_magma_gene_z.sh` to point at your binary.

### Step 2 — instruments, MR, meta-analysis, colocalisation (C1, C2, D1, D2, F1, G1)

```bash
cd ../02_instruments
python3 C1_build_instruments.py
python3 C2_harmonise.py

cd ../03_mr_meta_coloc
python3 D1-D2-F1-G1_mr_meta_coloc.py
```

Checkpoints: 1,981 instruments / 1,868 *cis*-eGenes → 1,371 genes / 1,430
instruments after harmonisation → 316 candidate *cis*-eGenes (MAGMA z ≥ 1), 225 with a
discovery estimate → 13 at BH-FDR < 0.05 → 9 direction-consistent at meta-P < 0.05 →
6 Tier 1, 7 retained after merging with the colocalisation tiers.

### Step 3 — overlap correction (F2)

```bash
cd ../04_overlap_correction
python3 F2_overlap_correction.py
```

Reads `results/figure4/fig4_mr_gene.csv` and writes
`results/figure4/fig4_overlap_correction.csv`. Expected output:

```
[overlap] candidates = 316
[overlap] inflation factor = 1.191
[overlap] meta-P < 0.05: 64 -> 59
[overlap] losing significance (5): APOBEC2, CPNE1, HLA-DRB5, NCK2, PM20D1
```

### Step 4 — single-cell QC, integration, annotation, scDRS (I1, J1)

```bash
cd ../06_singlecell
python3 I1a_sc_qc_build_adata.py gse173193
python3 I1b_sc_integrate_harmony.py norm
python3 I1b_sc_integrate_harmony.py embed
python3 I1c_sc_annotate.py gse173193

cd ../07_scdrs
python3 J1a_scdrs_input.py
python3 J1b_scdrs_score.py
python3 J1c_scdrs_downstream.py gse173193
```

QC: ≥ 200 genes per cell, mitochondrial reads ≤ 20 %, genes detected in ≥ 3 cells;
3000 HVGs with `batch_key`; PCA 30 components; Harmony on `batch`; Leiden at
resolution 0.8; nine annotated compartments.

### Step 5 — MR-PheWAS (M1)

```bash
cd ../09_phewas
python3 M1_mr_phewas.py
```

Queries the FinnGen R12 Results API for every endpoint in the 2,469-endpoint manifest,
for the two Tier 1 lead variants. Expected: 4,938 gene–endpoint tests (2 genes × 2,469).

### Step 6 — sensitivity diagnostics and tables (E1, G2)

```bash
cd ../05_sensitivity
python3 E1_mr_sensitivity.py

cd ../10_tables
python3 G2_tier_flags_and_tables.py
```

`E1` writes `suppl5_multi_iv_diag.csv`, `suppl5_per_iv.csv` and `suppl5_meta.json`
(Additional file 10 / Figure S5). `G2` rebuilds supplementary tables 1–5 and asserts the
stage counts (316 / 13 / 9 / 7) before writing; see §6 for the inputs it expects from
the exploratory analyses.

---

## 6. What is not included

* **P1 exploratory analyses** — tissue specificity against whole-blood eQTLGen
  (`fig4_blood_mr.csv`, `fig4_panelC_matrix.csv`), cell-type expression panel,
  pseudobulk violin, Hallmark/TRRUST activity, GO/KEGG, and CellChat. These produce
  Additional files 5 (Table 5a/5c), 6–8 and Figure 6.
* **P2 scripts** — figure rendering and table layout.
* **Small intermediates consumed by `G2`** — `suppl1_qc_meta.json`,
  `gse282038.runA_condition.csv`, `suppl2_replication.csv`. Table 4 and Table 5b are
  self-contained; Table 1 and Table 5a/5c require files from the steps above.
* Raw sequencing data, `.h5ad` objects, VCF/PLINK binaries and intermediate caches.

Each excluded step is named above, so a reader can tell which reported numbers are
reproducible from this repository alone and which depend on the exploratory analyses.

---

## 7. Reproducibility notes

* **Instruments** — one lead variant per SuSiE credible set (highest PIP), `P < 1e-5`,
  SNV only, `F = (beta/SE)^2 > 10`, no additional LD pruning.
* **Harmonisation** — effect alleles aligned to the exposure; palindromic variants with
  MAF 0.40–0.60 removed; unresolvable alleles dropped.
* **MR** — Wald ratio for single-instrument genes, fixed-effect IVW otherwise;
  DerSimonian–Laird random effects across the two GWAS, with `tau^2` retained for the
  overlap correction.
* **Colocalisation** — Wakefield approximate Bayes factors (`coloc.abf` independent
  approximation, priors p1 = p2 = 1e-4, p12 = 1e-5) computed over the variants of each
  gene's SuSiE credible set; Tier 1 at PP.H4 > 0.80, Tier 2 at 0.50 < PP.H4 ≤ 0.80.
* **scDRS** — a single scoring pass over all 55,414 cells; compartment association is
  then summarised in all cells (Run A), control cells (Run B) and preeclampsia cells
  (Run C), with Benjamini–Hochberg recomputed within each subset rather than the score
  being re-derived three times.
* **Covariates** — library size (`log_n_genes`) and mitochondrial fraction (`pct_mt`)
  enter the scDRS model as covariates; the embedding and clustering use Harmony over
  donor without regressing those covariates out first.
* **Overlap correction** — replication SE × 1.19, `tau^2` fixed at the primary
  estimate, meta-P re-derived; 64 → 59 of 316 candidates.
* **PheWAS** — BH-FDR < 0.05 within each gene, Bonferroni P < 2.0e-5 reported
  alongside; per-1-SD-expression Wald estimates and per-effect-allele associations are
  kept as separate quantities and never mixed.

---

## 8. License

MIT — see `LICENSE`. Please cite the manuscript and this repository together.
