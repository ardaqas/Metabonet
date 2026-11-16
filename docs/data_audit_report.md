# Phase 0 Data Audit Report

**Project:** MetaboNet: predicting CCLE metabolite abundance from Human-GEM network-gene expression
**Date:** 2026-11-16 · **Scope:** Phase 0 only (`research_plan.md` §4). No model was fitted and no metabolite→gene mapping was done.
**Reproduce:** `make data && make test && make audit` (see `docs/data_sources.md`). Notebook: `notebooks/00_data_audit.ipynb`.

## 1. Verdict

| Go/no-go criterion (research_plan §4) | Result | Status |
|---|---|---|
| ≥ 400 solid-tumour lines with metabolomics, RNA-seq, lineage and medium | **735** | **PASS** |
| ≥ 90% medium assignment coverage among matched solid lines | **100%** (735/735) | **PASS** |
| ≥ 15 polar metabolites pass target selection | not evaluated | Phase 1 |

**Recommendation: proceed to Phase 1.** Every required source was obtained and matched, both Phase 0 gates pass with a wide margin, and none of the open issues blocks the analysis. Section 7 lists the issues to carry into Phase 1.

## 2. Data obtained

The full provenance is in `docs/data_sources.md` and `data/raw/MANIFEST.tsv` (SHA-256, MD5, HTTP Last-Modified/ETag).

| Source | Status |
|---|---|
| `CCLE_metabolomics_20190502.csv` | ✅ Broad CCLE server (2019-05-02) |
| CCLE 2019 RSEM gene TPM (`CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz`) | ✅ Broad CCLE server |
| `Cell_lines_annotations_20181226.txt` | ✅ Broad CCLE server |
| DepMap `Model.csv` | ✅ **DepMap Public 26Q3**, downloaded manually with the release `README.txt`. The README header names the release, and its documented Model.csv schema matches the file (48 of 48 columns, same order). The DepMap portal blocks automated access with a Cloudflare Turnstile check (decision D10). |
| Li et al. 2019 Supplementary Table 1 | ✅ Nature Medicine ESM (`41591_2019_404_MOESM2_ESM.xlsx`) |
| Cellosaurus (problematic lines, same-donor links) | ✅ API, release 57.0, 927 RRIDs, all resolved |

No source was substituted. The DepMap portal blocker is described in `data_sources.md`.

## 3. Cell-line matching

| Step | Lines |
|---|---:|
| CCLE 2019 metabolomics lines | 928 |
| … with a DepMap ID | 928 |
| … found in current Model.csv (ModelID or alias) | 928 |
| … with CCLE 2019 RNA-seq | 913 |
| … with lineage (OncotreeLineage) | 913 |
| … duplicate cell lines removed | 912 |
| … cancer models (16 non-cancerous lines removed) | 896 |
| … solid tumour (161 Lymphoid/Myeloid lines removed) | 735 |
| … with medium assigned = **primary set** | **735** |

- **RNA-seq overlap.** 913 of 928 metabolomics lines have CCLE 2019 RNA-seq. The other 15 have no column under any name or ID: BECKER, CGTHW1, CHL1, CMK115, COLO704, DOV13, MOTN1, NCIH2141, OC315, OC316, OUMS27, OVCAR5, SHSY5Y, SKMEL2, TK10.
- **Spot check.** 20 random primary-set lines (seed 0) were checked by hand: CCLE_ID ↔ DepMap ID ↔ Model.csv CCLEName ↔ Cellosaurus name ↔ lineage ↔ medium. All were consistent.
- The **frozen table** is `data/processed/cell_lines_phase0.csv`: 928 rows, 44 columns, with the `in_primary_set` flag. Hashes are in `data/processed/FREEZE.tsv`.

## 4. Identifier and cell-line issues

All issues are listed one per line in `data/processed/id_issues.csv`.

| Issue | n | Handling |
|---|---:|---|
| Missing DepMap_ID in the metabolomics file (OC315_OVARY) | 1 | Filled from the 2019 annotation (ACH-001144). The line has no RNA-seq anyway. |
| DepMap ID now a retired alias: **HCC1588_LUNG** (ACH-001078 → ACH-000007 = **LS513_LARGE_INTESTINE**, which is also in the set) | 1 | **Duplicate cell line.** HCC1588 excluded, LS513 kept (D16). |
| Name changed in current Model.csv | 4 | KE97 → `KE97_STOMACH` (Lymphoid, excluded as haematopoietic). **NCIH1339 → NCI-H157-DM** (ICLAC: the DepMap stock is NCI-H157; still Lung, kept, flagged). NCIH684_LIVER → `_LARGE_INTESTINE` (lineage Bowel used). HCC1588 (above). |
| Li et al. names differ from CCLE 2019 IDs | 15 | 14 lines were renamed `*_FIBROBLAST` after Li et al. 2019, and S117 moved from SOFT_TISSUE to THYROID. Matched by unique cell-line name. |
| Non-cancerous models (14 fibroblasts renamed after Li et al., HS729 and EHEB) | 16 | Excluded from the solid set (D11). These 14 fibroblast lines were sold as tumour lines but are now known to be fibroblasts. |
| Cellosaurus problematic lines | 51 (**26 in primary set**) | Kept and flagged (D15). In the primary set: 8 misclassified, 8 contaminated, 3 partially contaminated and 7 possibly/probably misclassified or misidentified. |
| Same-donor groups (DepMap PatientID, RRID, Cellosaurus OI/HI/derivative) | 20 groups, 41 lines (**11 groups, 22 lines in primary set**) | `donor_group_id` must be one CV group (D14). Examples: SW480/SW620, KPL1/MCF7, ONCODG1/OVCAR3, PATU8988S/T, IGR37/IGR39, HEC1A/B. |
| Haematopoietic flag disagrees with 2019 Site_Primary | 1 | TO175T (2019: haematopoietic; current: Fibroblast). Current lineage used; non-cancerous anyway. |
| Duplicate CCLE_ID or DepMap ID in the metabolomics or RNA-seq matrices | 0 | — |

## 5. Lineage and growth type (primary set, n = 735)

- **Lineages, raw OncotreeLineage (22):** Lung 169, Esophagus/Stomach 63, Bowel 56, Skin 49, CNS/Brain 49, Breast 48, Ovary/Fallopian Tube 45, Pancreas 41, Head and Neck 33, Uterus 30, Kidney 27, Bladder/Urinary Tract 25, Liver 24, Peripheral Nervous System 16, Bone 15; then Soft Tissue 10, Thyroid 10, Pleura 9, Prostate 7, Biliary Tract 6, Ampulla of Vater 2, Cervix 1.
- **After collapsing groups of fewer than 15 lines:** 16 categories, with **Other = 45**.
- **Growth pattern (Model.csv):** Adherent 633, Mixed 55, Suspension 44, Unknown 3. All are kept and flagged (D11). **This is not the population of Li et al. 2022**, who used 454 solid *adherent* lines. Any comparison with their numbers must say so (D21).

## 6. Culture medium

- **Coverage:** Li et al. Supp. Table 1 gives a medium for all 928 lines (42 distinct strings), so **100%** of the primary set is assigned from it. The 2019 `Growth.Medium` alone would cover 74.3% (546/735). It is blank for 224 of the 928 metabolomics lines and has 216 different spellings.
- **Agreement:** where both sources can be classified (546 lines), they agree for **99.5%** (543). The three disagreements are:
  - 253J: Li "DMEM", 2019 "Eagle";
  - GB1: Li "EMEM", 2019 "Dulbecco's modified Eagle's medium";
  - NCIH684: Li "DMEM", 2019 "ACL-4".

  Li is used in all three cases (D12).
- **Base media in the primary set:**
  - RPMI-1640: 375
  - DMEM: 114
  - EMEM/MEM: 95
  - DMEM:F12: 59
  - McCoy's 5A: 26
  - mixtures: 14 (RPMI+MEM 4, RPMI+F12 4, MCDB105+M199 3, EMEM+F12 2, DMEM+RPMI 1)
  - Leibovitz L-15: 12
  - IMDM: 11
  - Ham's F12: 9
  - Ham's F10: 6
  - Waymouth's: 5
  - Ham's F12K: 4
  - Williams' E: 4
  - alpha-MEM: 1
- **Model categories (fewer than 15 lines → other):** RPMI-1640 375, DMEM 114, EMEM/MEM 95, **other 66**, DMEM:F12 59, McCoy's 5A 26.
- **Normalisation table:** `data/processed/medium_normalization.csv` has one row for each unique raw string and source, giving its rule, category and line counts. Every raw string was reviewed by hand. The two-source cross-tabulation is in `data/processed/medium_source_agreement.csv`. The rules have 49 unit tests (`tests/test_medium.py`).
- **Caveat:** medium is collinear with lineage. As planned, the two are modelled jointly and their separate coefficients are not interpreted.

## 7. Metabolomics

- **Matrix:** 928 lines × 225 metabolites. The released CSV is **identical** to Li et al.'s `1-clean data` sheet (maximum absolute difference 0).
- **Transformation:** values are on a **log10 scale**.
  - Per-metabolite SD is 0.36 (median), close to log10(raw) at 0.41 and far from log2 (1.37) or ln (0.95).
  - Values correlate better with log10(raw peak area) (median r = 0.90) than with linear raw areas (r = 0.73).
  - Every metabolite is re-centred near about 5.9 (overall range 2.98–9.10). Values can therefore be compared across lines **within** a metabolite, not between metabolites.
  - Median normalisation and batch correction give a regression slope on log10(raw) of 0.77 (IQR 0.68–0.86).
  - **Do not log-transform again.**
- **Missing values:** the released CSV has none, because the data generators imputed with MICE. Missingness before imputation, recovered from `1-raw data`, is 19,185 of 208,800 cells (**9.2%**), matching the paper's "~9%". All of those cells are imputed in the CSV.
  - Polar metabolites in the primary set: median 2.0% missing, 48 with none, **21 above 20%**, 8 above 30%. The highest are α-hydroxybutyrate (37%), acetylcholine (36%), taurocholate (36%), sarcosine, cAMP and hexoses (HILIC pos) (≈34–35%).
  - Per-metabolite values are in `metabolite_inventory.csv`.
  - **Rule fixed (D18):** a metabolite is eligible only if at most 20% of its primary-set values were missing before imputation. **115 of 136 polar metabolites pass**, including all five positive controls. The flag is `missingness_rule_pass`.
- **Technical CV (Li):** 221 of 225 metabolites have CV < 0.1. The four above 0.1 are all polar: acetylcholine 0.21, pipecolic acid 0.19, serotonin 0.14, lactose 0.13.
- **Assay class (`metabolite_inventory.csv`):**
  - From Li's LC-MS parameter sheets: **136 polar (HILIC; 77 positive and 59 negative ion mode) + 89 lipid (RPLC)**.
  - 8 polar data columns use short names and are matched to the sheet through an explicit alias table. `DHAP/glyceraldehyde 3P` is one data feature with two transitions, which accounts for the sheet's 137 rows.
  - **Discrepancy:** the paper's text says 124 polar + 101 lipid. The method-sheet counts reproduce Cavicchioli et al.'s 136/89. The 14 carnitines/acylcarnitines are measured by HILIC, which may explain most of the gap if the paper counted them as lipids, but this is **not** confirmed.
  - As the plan requires, the assay method sheet is used.
- **Features named as possibly more than one species:** 13 polar features, e.g. `F1P/F6P/G1P/G6P`, `succinate/methylmalonate`, `SDMA/ADMA`, `cis/trans-hydroxyproline`. They are flagged here only; Phase 1 mapping decides about them.

## 8. Expression (CCLE 2019 RSEM gene TPM)

- **Shape:** 57,820 genes × 1,019 lines, with no duplicate line columns. 913 of the 928 metabolomics lines are present.
- **Identifiers:** versioned GENCODE v19 Ensembl IDs. 57,773 are `ENSG…` and **47 are `ENSGR…`**, the chrY PAR copies. No IDs collide after the version suffix is stripped. All 47 `ENSGR` rows have an `ENSG` chrX counterpart and are zero in every line. **They are dropped (D22)**, leaving 57,773 genes. The file also has a `transcript_ids` column, which is ignored.
- **Values:** **linear TPM, not log**. Column sums are 999,995–1,000,006. There are no NaNs and no negative values, 54% of entries are zero, and the maximum is 225,730.
- **Expression filter (for information only):** 15,106 genes have mean log2(TPM+1) ≥ 1 in the primary set, and 8,769 are zero in every primary line. The protein-coding restriction and the final filter are set in Phase 1/2. The GENCODE v19 biotype annotation was not downloaded, because Phase 0 doesn't need it.

## 9. Unresolved items for Phase 1

1. **Polar/lipid count discrepancy** (136/89 vs 124/101). It doesn't affect Phase 1, which classifies by method sheet as planned. Note it in the write-up.
2. ~~Pre-imputation missingness threshold~~ **Resolved (D18):** at most 20% missing, written into `expectations.md` §1a. 115 of 136 polar metabolites remain eligible.
3. ~~Batch labels~~ **Resolved as a limitation (D19):** no source has them. The batch sensitivity analysis is dropped and logged as a deviation in `expectations.md` §8.
4. **Gene ID mapping:** GENCODE v19 → current Human-GEM Ensembl IDs (Phase 1). The `ENSGR` rows are already settled: they are dropped (D22).
5. ~~Problematic lines~~ **Resolved (D20):** a fourth sensitivity analysis excludes the 26 flagged lines (n = 709). It was added before any results and is logged as a deviation.
6. **DepMap portal byte-identity:** the CCLE 2019 files from the Broad server could not be compared with the portal copies.
7. ~~Model.csv release~~ **Resolved.** The release `README.txt` header says DepMap Public 26Q3, its documented Model.csv schema matches the file exactly, and `make verify` now checks this on every run. One limit remains: the README has no checksum for Model.csv, so this confirms the release from the schema, not from the bytes.

## 10. Outputs

| File | Content |
|---|---|
| `data/processed/cell_lines_phase0.csv` | Frozen matched cell-line table: IDs, lineage (raw and collapsed), flags, medium, donor groups, `in_primary_set` |
| `data/processed/medium_normalization.csv` | One row per unique raw medium string and source: category, rule, counts |
| `data/processed/medium_source_agreement.csv` | Li vs 2019 annotation cross-tabulation |
| `data/processed/metabolite_inventory.csv` | Per metabolite: assay class, ion mode, CV, missingness, spread |
| `data/processed/id_issues.csv` | Every matching issue and how it was handled |
| `data/processed/attrition.csv` | Line counts at each step |
| `data/processed/audit_summary.json` | All statistics in this report |
| `data/processed/FREEZE.tsv` | SHA-256 of the files above (identical across repeated runs) |
| `docs/data_sources.md`, `docs/decision_log.md` | Provenance; decisions D10–D17 |
