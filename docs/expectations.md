# Expectations: Testable Hypotheses and Evaluation Criteria

**Project:** Predicting metabolite abundance from gene expression using human metabolic network information (CCLE)
**Status:** Draft v0.1 (16 November 2026). **To be frozen** (dated, committed, with the final target list appended) at the end of Phase 1 and before any model is fitted. After freezing, changes are allowed only as logged deviations (Section 8).

This document states in advance what we expect, how each expectation will be tested, and what each outcome would mean. Its purpose is to prevent post-hoc reinterpretation. It does not assume the main hypothesis is true; based on the literature, a mixed or largely negative result is the most likely outcome.

---

## 1. Definitions

- **Target:** a polar CCLE metabolite that passes the outcome-blind selection rules in `research_plan.md`, Section 5.
- **Gₘ:** first-neighbour network gene set of target *m* (size *k*), from the pinned Human-GEM release.
- **M0:** lineage + medium baseline. **M1:** M0 + Gₘ. **M2:** M0 + all expressed genes.
- **N-met / N-all:** M0 + *k* random genes from Human-GEM metabolic genes / from all expressed genes (B = 1,000 draws each).
- **P:** M1 fitted to the target permuted within lineage × medium strata.
- **R²:** out-of-fold R² (5-fold × 10 repeats, identical folds across models), averaged over repeats.
- **ΔR²(X) = R²(M0 + X) − R²(M0).**
- **Empirical p (vs N-met or N-all):** (1 + #{random ΔR² ≥ observed ΔR²(M1)}) / (B + 1).
- **Minimum meaningful effect:** ΔR² ≥ 0.02 (2 percentage points of variance). Chosen because smaller gains are within typical repeat-to-repeat variation for n ≈ 500 and would not change biological interpretation. Results are reported for all targets regardless of this threshold.

### 1a. Pre-specified data rules (fixed 2026-11-16, after the Phase 0 audit and before any mapping or modelling)

- **Primary cell-line set:** the 735 solid-tumour lines with `in_primary_set = True` in `data/processed/cell_lines_phase0.csv` (frozen; hashes in `data/processed/FREEZE.tsv`). Definition: `docs/decision_log.md` D11 and D16.
- **Missingness rule (target rule 4):** a metabolite is eligible only if **≤ 20 %** of its values in the primary set were missing **before** the data generators' MICE imputation (from Li et al. 2019 Supp. Table 1, `1-raw data`). Above 20 %, a fifth or more of the values are imputed guesses that borrow information from other lines. The rule removes 21 of 136 polar metabolites, leaving 115 eligible on this criterion. All five literature positive controls pass (maximum: kynurenine, 11 %). Per-metabolite flag: `missingness_rule_pass` in `data/processed/metabolite_inventory.csv`.
- **Expression genes:** the 47 GENCODE v19 `ENSGR…` rows (chrY pseudoautosomal copies) are dropped before any gene filtering or ID mapping. Each has an `ENSG…` chrX counterpart, and every one is zero in all 1,019 lines.
- **Lineage and medium categories:** `lineage_collapsed` and `medium_category` from the frozen table (categories with fewer than 15 primary-set lines become "Other"/"other").
- **CV grouping:** `donor_group_id` from the frozen table. Lines in the same donor group are always placed in the same fold.

---

## 2. Sanity checks that must pass before any hypothesis is evaluated

| ID | Check | Pass criterion | If it fails |
|---|---|---|---|
| S1 | Permuted-target control (P) | Median ΔR² across targets within ±0.005 and no target with ΔR² > 0.02 | Leak or bug; stop and debug |
| S2 | Synthetic target recovery | A simulated target built as a noisy linear function of 5 genes in a defined set is recovered by M1 with ΔR² close to the simulated signal, and beats N-met for ≥ 95% of simulations | Pipeline cannot detect true signal; fix before proceeding |
| S3 | Leakage checklist (`research_plan.md`, Section 10.3) | All items "no" | Fix; rerun |
| S4 | Data audit | Go/no-go criteria in `research_plan.md`, Section 4 met | Narrow the question (e.g., single medium) rather than add complexity |

---

## 3. Hypotheses

### H1 — Network genes add information beyond lineage and medium

- **Rationale:** Enzymes and transporters that produce or consume a metabolite are its most direct transcriptional determinants (Zelezniak et al. 2014); several single-enzyme links are established in CCLE.
- **Prediction:** ΔR²(M1) > 0 for most targets.
- **Test:** One-sided Wilcoxon signed-rank test of ΔR²(M1) across targets against 0; count of targets with ΔR²(M1) ≥ 0.02.
- **Supported if:** Wilcoxon p < 0.05 **and** at least 25% of targets have ΔR²(M1) ≥ 0.02.
- **Interpretation caution:** Support for H1 alone says nothing about network specificity — any gene set may add information (see H2).

### H2 (primary) — Network genes beat size-matched random metabolic genes

- **Rationale:** If metabolite levels are partly controlled locally, Gₘ should be more informative than arbitrary metabolic genes. Competing evidence: random genes performed about as well as metabolic genes in Recon8D, and proximal pairs were only a small fraction of gene–metabolite associations in tumours (Benedetti et al. 2023).
- **Prediction:** For a minority of targets, ΔR²(M1) exceeds the N-met null; across targets, ΔR²(M1) exceeds the median N-met ΔR² more often than expected by chance.
- **Tests:**
  - Global: one-sided Wilcoxon signed-rank test of ΔR²(M1) − median ΔR²(N-met) across targets (α = 0.05); binomial test of the number of targets above their own 95th N-met percentile vs 5% expected.
  - Per target: empirical p vs N-met, Benjamini–Hochberg FDR q < 0.10, and ΔR²(M1) ≥ 0.02. Targets meeting both are called **network-specific**.
- **Supported if:** global Wilcoxon p < 0.05 **and** at least 3 targets are network-specific.
- **Partially supported if:** 1–2 network-specific targets (typically positive controls) but no global effect. Interpretation: local transcriptional control is detectable for specific metabolites, not a general rule.
- **Not supported if:** no network-specific targets and global p ≥ 0.05. Interpretation: in CCLE, expression-based metabolite prediction mostly reflects general cell state.

### H3 — Network genes also beat random genes from the whole transcriptome

- **Test:** As H2, against N-all.
- **Expected relation to H2:** Beating N-all is easier than beating N-met, because metabolic genes as a class carry tissue and cell-state information (Peng et al. 2018). If M1 beats N-all but not N-met, the conclusion is "metabolic genes are informative", not "network locality is informative".

### H4 — The all-gene model is at least as good as the network model

- **Rationale:** Li et al. (2022) found all-gene models generally performed better than metabolic subsets; flexible all-gene models capture lineage-correlated cell state.
- **Prediction:** ΔR²(M2) ≥ ΔR²(M1) for most targets.
- **Test:** Paired Wilcoxon signed-rank of ΔR²(M2) − ΔR²(M1) across targets; capture ratio ΔR²(M1)/ΔR²(M2) for targets with ΔR²(M2) > 0.02.
- **Expected result:** M2 better on average. The interesting outcome is the subset of targets where the capture ratio is high (e.g., ≥ 0.7), i.e., a small, interpretable gene set recovers most of the achievable signal.
- **Surprising (and to be scrutinised for leakage) if:** M1 substantially beats M2 for many targets.

### H5 — Literature positive controls show network-specific signal

- **Targets (if they pass selection rules):** 1-methylnicotinamide (NNMT), kynurenine (IDO1/TDO2), creatine (SLC6A8/GATM), N-carbamoyl-β-alanine (DPYD), taurine (SLC6A6).
- **Prediction:** Each control that passes selection has ΔR²(M1) ≥ 0.02 and empirical p vs N-met < 0.05 (uncorrected; these are confirmatory).
- **If controls fail:** Either the pipeline is too conservative (check S2) or the two-stage lineage/medium adjustment absorbs the signal (e.g., MNA is lineage-associated). Report both M0-adjusted and unadjusted results for controls to diagnose which.

### H6 — Target class matters

- **Prediction:** Medium-supplied nutrients (e.g., essential amino acids) show higher R²(M0) relative to ΔR²(M1) and fewer network-specific calls than synthesised intermediates.
- **Test:** Compare ΔR²(M1) − median ΔR²(N-met) between pre-annotated classes (Mann–Whitney U); descriptive given small group sizes.
- **Status:** Exploratory. Reported regardless of direction.

---

## 4. Honest prior expectations

Stated so they can be compared with outcomes, not to steer the analysis.

- R²(M0) will be substantial for many targets (lineage effects reported for most metabolites by Li et al. 2019).
- ΔR²(M1) will be small for most targets (median plausibly below 0.05).
- Random gene sets will **not** have ΔR² ≈ 0; they will show modest positive gains because expression encodes cell state. The null is not zero.
- Network-specific signal is expected mainly for metabolites with a dedicated, transcriptionally variable enzyme or transporter (positive-control type). A plausible range is 10–25% of targets; fewer would not be surprising.
- M2 will outperform M1 on average.

---

## 5. Interpretation matrix

| H1 | H2 (vs N-met) | H3 (vs N-all) | Conclusion |
|---|---|---|---|
| No | No | No | Beyond lineage and medium, expression carries little held-out information about these metabolites in CCLE |
| Yes | No | No | Expression helps, but network genes are no better than random genes: signal is general cell state |
| Yes | No | Yes | Metabolic genes are informative as a class, but local network position adds nothing detectable |
| Yes | Partial | Yes | Local network control is detectable for specific metabolites (likely those with dedicated enzymes), not as a general rule |
| Yes | Yes | Yes | Network locality carries specific predictive information across many metabolites (would contradict the random-gene result of Recon8D; check carefully for leakage and target-selection bias) |

Every row is a publishable-quality, honest outcome for a portfolio project. The project succeeds if the design can distinguish these rows, not if a particular row is obtained.

---

## 6. Sensitivity analyses (pre-specified; maximum four)

1. All cell lines including haematopoietic/suspension lines.
2. Enzymatic-only gene sets (transport reactions excluded).
3. Variance-matched random gene sets.
4. Primary set **excluding the 26 lines that Cellosaurus flags as problematic** (misclassified, contaminated or misidentified; `is_problematic` in the frozen table; n = 709). Added 2026-11-16, before any mapping or modelling; see the deviation log.

A conclusion is called **robust** if the H2 classification (supported / partial / not supported) is unchanged in all four.

The batch-covariate sensitivity analysis considered in `research_plan.md` §2.3 is **not possible**: no public source gives the LC-MS batch labels (Phase 0 audit). It is reported as a limitation and is not replaced.

---

## 7. Reporting commitments

- Report results for **every** target, including failures and positive controls that fail.
- Report R²(M0), ΔR²(M1), ΔR²(M2), the N-met and N-all null distributions, empirical p-values and FDR q-values per target.
- Report the number of cell lines, metabolites and genes lost at every matching and filtering step.
- Report the metabolite mapping table with exclusion reasons.
- State limitations explicitly: different cultures for RNA and metabolites; relative abundances; medium heterogeneity; mRNA as a proxy for enzyme activity; cell lines vs tumours; associative, not causal, inference; **no LC-MS batch labels are available, so residual batch structure across the four batches cannot be checked or adjusted**.
- When comparing with Li et al. (2022), state that their 454 lines were solid **and adherent**. The primary set here is solid tumours of any growth type: 633 adherent, 55 mixed, 44 suspension and 3 unknown. The two populations are therefore not identical.
- Do not describe results as "predicting metabolism" in general; describe them as held-out prediction of relative abundance of the specific targets studied.

---

## 8. Deviation log

Any change after freezing (targets, gene-set rules, thresholds, model, CV design) is recorded here with date, change, reason, and whether it was made before or after seeing results.

| Date | Change | Reason | Before/after seeing results |
|---|---|---|---|
| 2026-11-16 | Target rule 4 ("missingness acceptable") made specific: ≤ 20 % missing before imputation in the primary set (§1a) | The plan required a threshold; the Phase 0 audit made missingness measurable | **Before** (pre-freeze, set before mapping; no model fitted) |
| 2026-11-16 | Batch-covariate sensitivity analysis (research_plan §2.3) dropped and recorded as a limitation | Batch labels are not in any public source examined (Phase 0 audit) | **Before** (no model fitted) |
| 2026-11-16 | Fourth sensitivity analysis added: exclude 26 Cellosaurus-problematic lines; the "maximum three" limit is raised to four | Expected reviewer question about misidentified or contaminated lines. The three-analysis limit was meant to stop choices made after seeing results, and no results exist yet | **Before** (no model fitted) |
| 2026-11-16 | 47 `ENSGR` chrY pseudoautosomal gene rows dropped from expression | Duplicates of chrX genes; all zero in the data | **Before** (no model fitted) |

---

## 9. Target list (to be appended at freeze)

To be filled at the end of Phase 1 with: target name, HMDB ID, Human-GEM metabolite ID (compartment-free), gene-set size *k*, class (synthesised / medium-supplied / mixed), positive-control flag. Frozen on: ______ (commit: ______).
