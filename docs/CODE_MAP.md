# Code map: script → manuscript → expected numbers

Every script is listed with the section of the manuscript it corresponds to, the
parameters that can be checked against the text, and the files it writes. Numbers in
**bold** are the checkpoints reported in the Results; a run is consistent with the
manuscript only if they are reproduced.

---

## 1. Instruments and harmonisation

| Script | Manuscript | Parameters to check | Output / checkpoint |
|---|---|---|---|
| `02_instruments/C1_build_instruments.py` | Methods, *Mendelian randomisation…* | lead variant per credible set = highest PIP; `P < 1e-5`; SNV only; `F = (beta/SE)^2 > 10`; no LD pruning | `fig4_instruments.csv`; **1,981 instruments / 1,868 *cis*-eGenes** |
| `02_instruments/C2_harmonise.py` | Methods, same section; Table 2 | effect allele aligned to exposure; palindromic variants with MAF 0.40–0.60 dropped; unresolvable alleles dropped | `fig4_outcome_finngen.csv`, `fig4_outcome_gwas2.csv`; **1,371 genes / 1,430 instruments** |

`C2` queries GCST90269903 through the EBI tabix index when the local harmonised file is
absent; the GB download is 375 MB, the tabix route fetches ±1 bp per instrument.

---

## 2. MR, meta-analysis, colocalisation

| Script | Manuscript | Parameters to check | Output / checkpoint |
|---|---|---|---|
| `03_mr_meta_coloc/D1-D2-F1-G1_mr_meta_coloc.py` | Methods, *Mendelian randomisation…*; Results, Fig. 4A | Wald ratio (1 instrument) / fixed-effect IVW (≥ 2); `OR = exp(beta)`; DerSimonian–Laird across GWAS; BH-FDR within the **225**-gene family | `fig4_mr_gene.csv`, `fig4_candidate_genes.csv`, `fig4_funnel.csv`; **316 → 13 → 9** |
| same file, colocalisation block | Methods, *Colocalisation…*; Table 3 | `coloc.abf` Wakefield ABF, priors p1 = p2 = 1e-4, p12 = 1e-5; variants = SuSiE credible set of the gene; Tier 1 > 0.80, Tier 2 0.50–0.80 | **Tier 1 = 6, Tier 2 = 9** among 316 candidates |
| `04_overlap_correction/F2_overlap_correction.py` | Methods, *overlap correction*; Results | replication SE × √(296,824/209,254) = 1.191; `tau^2` fixed at the primary estimate | **meta-P < 0.05: 64 → 59**; five genes lost: APOBEC2, CPNE1, HLA-DRB5, NCK2, PM20D1; largest inflated meta-P among the nine = **0.037** (ELP3) |
| `05_sensitivity/E1_mr_sensitivity.py` | Results, Additional file 10 / Fig. S5 | Cochran's Q; MR-Egger via weighted regression of `b_out` on `b_exp`; leave-one-out; Steiger `r2 = z^2/(z^2+N)` | `suppl5_multi_iv_diag.csv`, `suppl5_per_iv.csv`; **1,430 instruments, 1,313 single-instrument genes, 58 testable loci, 2 with heterogeneity** |

Colocalisation is computed **over the variants of each gene's credible set**, not over a
fixed megabase window: `coloc_n_SNP` in Additional file 4 ranges from 1 to 363 and is 3
for ERAP2 and 19 for TK2. The manuscript's Methods sentence states this directly (median
10 variants per gene, range 1–363), so the text, the table and this script describe the
same variant set.

---

## 3. MAGMA and compartment sets

| Script | Manuscript | Parameters to check | Output |
|---|---|---|---|
| `01_magma/K1a_magma_input_pval.py` | Methods, *Gene set enrichment* | FinnGen R12, N = 268,336 | `finngen_O15_PREECLAMPS.magma.pval` |
| `01_magma/K1b_magma_gene_z.sh` | same | MAGMA v1.10; 1000G EUR; `--annotate window=10,10`; NCBI37.3 | `finngen_O15.genes.out` |
| `01_magma/K1c_magma_genesout_to_pval.py` | same; Methods, *scDRS* | entrez → symbol; top-1,000 gene set for scDRS | `finngen_O15.pval.tsv` |
| `08_compartments/J2_compartment_magma_sets.py` | Methods, *Compartment…* | the nine marker sets are the same curated lists used for annotation in `I1c` | `magma_celltype_sets.txt` |

The gene-level z used to build the candidate pool is the same quantity in two places:
`z ≥ 1` gives the **316**-gene pool (225 with a discovery estimate) and `z ≥ 3` gives
the **179**-gene universe used for pathway analysis.

---

## 4. Single cell

| Script | Manuscript | Parameters to check | Output / checkpoint |
|---|---|---|---|
| `06_singlecell/I1a_sc_qc_build_adata.py` | Methods, *Placental single-cell resource* | ≥ 200 genes/cell; mitochondrial ≤ 20 %; gene in ≥ 3 cells (per sample) | `combined_raw.h5ad`; **55,414 cells from 58,577** |
| `06_singlecell/I1b_sc_integrate_harmony.py` | same | 3000 HVGs (`batch_key`), PCA 30, Harmony on `batch`, Leiden 0.8, UMAP | `combined_lognorm/hvg/emb.h5ad` |
| `06_singlecell/I1c_sc_annotate.py` | same | `sc.tl.score_genes` on nine canonical marker sets; cluster label = arg-max of the mean score | `combined_annotated.h5ad`; **9 compartments** |
| `07_scdrs/J1a_scdrs_input.py` | Methods, *Compartment disease relevance scoring* | raw counts kept for scDRS | `combined_scdrs.h5ad` |
| `07_scdrs/J1b_scdrs_score.py` | same | MAGMA top-1,000 set; 1,000 matched control sets; `adj_prop = cell_type`; covariates log n_genes and %mt | `O15_PREECLAMPS.full_score.gz` |
| `07_scdrs/J1c_scdrs_downstream.py` | same; Fig. 2B, C, E | compartment Monte-Carlo association (q95 vs control q95); BH within each subset; donor aggregation | `results/gse173193/runA_*, runB_*, runC_*, donor.csv`; **12,046 control cells, 15,846 PE cells** |

`J1c` **subsets** the cell-level scores from a single scoring run (`J1b`) — it does not
re-score the control-only or case-only cells. "Run B" and "Run C" therefore denote
compartment summaries computed within those cell subsets, with their own BH family, not
three independent scDRS runs.

---

## 5. Phenome-wide scan

| Script | Manuscript | Parameters to check | Output / checkpoint |
|---|---|---|---|
| `09_phewas/M1_mr_phewas.py` | Methods, *Phenome-wide association scan*; Fig. 5C | FinnGen R12 Results API; 2,469 endpoints; Wald ratio; SE recovered from `beta` and `pval` | `fig5C_phewas.csv`; **4,938 tests = 2 genes × 2,469**, BH within gene |

Only the two Tier 1 lead variants are screened (`ERAP2` rs4869313 `5:96888176:T:G`,
`TK2` rs12448207 `16:66598482:T:C`). The script asserts that the API alternate allele
equals the exposure effect allele before computing the ratio.

---

## 6. Tables

| Script | Manuscript | Notes |
|---|---|---|
| `10_tables/G2_tier_flags_and_tables.py` | Additional files 1–5 | Rebuilds tables 1–5 from the analysis products and asserts the stage counts (316 / 13 / 9 / 7) before writing; aborts if any assert fails. |

Inputs still required from the exploratory analyses (not in this repository):
`fig4_blood_mr.csv`, `fig4_panelC_matrix.csv`, `suppl1_qc_meta.json`,
`gse282038.runA_condition.csv`, `suppl2_replication.csv`. Table 3 and Table 4 need only
`fig4_mr_gene.csv` and `fig4_candidate_genes.csv`.

---

## 7. Stage counts

| Stage | n | Where |
|---|---|---|
| *cis*-eGenes with a lead variant (instruments) | 1,868 | `C1` |
| genes / instruments after harmonisation | 1,371 / 1,430 | `C2` |
| candidate pool (MAGMA z ≥ 1) | 316 | `D1` |
| with a discovery GWAS estimate (FDR family) | 225 | `D1` |
| discovery BH-FDR < 0.05 | 13 | `D1` |
| direction-consistent and meta-P < 0.05 | 9 | `D1` |
| colocalisation Tier 1 / Tier 2 | 6 / 9 | `D1` |
| retained after intersecting with the cascade | 7 | `G2` |
| cells after QC | 55,414 | `I1a` |
| PheWAS tests | 4,938 | `M1` |
