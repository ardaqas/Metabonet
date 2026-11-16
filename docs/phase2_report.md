# Phase 2 Report: Pipeline and Controls

**Project:** MetaboNet · **Date:** 2026-11-16 · **Scope:** Phase 2 only (`research_plan.md` §7–§9, §13; `expectations.md` §2).
**No hypothesis was evaluated.** The main analysis (Phase 3) has not been run. The only real-target fit in this phase is the pre-designated timing target (S5), whose result arrays are stored and were not used for any decision.
**Reproduce:** `make phase2` (folds are frozen in `results/phase2/folds.csv` and are not regenerated if present). Unit tests: `tests/test_phase2.py` (11 tests; 82 in the whole suite).

---

## 1. Verdict

| Check (expectations §2) | Criterion | Result | Status |
|---|---|---|---|
| S1 Permuted-target control | median ΔR²(P) across targets within ±0.005; no target > 0.02 | median **−0.0004**; max **0.008**; 0 of 74 above 0.02 | **PASS** |
| S2 Synthetic-target recovery | ΔR²(M1) close to the injected signal (median abs. error ≤ 0.05); beats N-met in ≥ 95% of simulations | median abs. error **0.045**; **20/20** beat N-met (all p = 0.005, the minimum at B = 200) | **PASS** |
| S3 Leakage checklist | all items "no" | all "no" with evidence (Section 5) | **PASS** |
| S4 Data audit | Phase 0 gates | passed in Phase 0 | PASS |
| S5 Compute contingency | projected total ≤ ~8 h at B = 1,000 | **1.7 h** projected (65 s per median-size target) | **PASS; B stays 1,000** |

**Recommendation: the pipeline is ready for Phase 3.** The one caveat worth carrying forward is in S2: the pipeline recovers about 70% of an injected linear signal, so ΔR² values in Phase 3 should be read as conservative estimates (Section 4).

## 2. What was built

`src/phase2/` implements research_plan §7–§9 exactly as written, with every learned quantity confined to the training fold:

| Component | Implementation | Decision |
|---|---|---|
| Inputs | 735 primary lines; 74 frozen targets + 1 diagnostic (MNA); log2(TPM+1) expression, ENSGR dropped; universes U_all = 12,666 genes (expressed protein-coding + expressed Human-GEM-mapped), U_met = 2,123 | D40 |
| Folds | `StratifiedGroupKFold`, 5 folds × 10 repeats, stratified by `lineage_collapsed`, grouped by `donor_group_id`, seeds 1000–1009; frozen in `results/phase2/folds.csv` | D42 |
| Baseline M0 | OLS on intercept + one-hot lineage + one-hot medium (22 columns), fitted per training fold; residuals carried to stage 2 | §7 |
| Gene features | standardised with training-fold mean/SD, computed once per fold and shared by all models | §7 |
| Ridge | closed-form, unpenalised intercept, penalty by exact LOO over 17 log-spaced values 10⁻²–10⁶ on the training fold; identical to `RidgeCV(cv=None)` (verified) | D41 |
| Prediction | test = M0 prediction + ridge prediction; R² pooled over the 5 folds of a repeat, averaged over 10 repeats; Spearman ρ secondary | §9 |
| Random sets | B draws of size k from U_met \ G_m (N-met) and U_all \ G_m (N-all), seeded per target, identical folds | §6 |
| P control | one permutation per target within lineage × medium strata, identical folds and procedure | D43 |
| Entry point | `run_target(name)` returns M0, M1, M1-enzymatic, M2, N-met, N-all and P arrays; Phase 3 calls it unchanged | — |

**Fold characteristics.** Fold sizes 145–149. No donor group is split (11 multi-line groups, 22 lines). The smallest lineage (Bone, 15 lines) is absent from 3 of the 50 test folds; in those folds all Bone lines are in training, so the baseline is still estimable and no test line lacks a lineage coefficient.

## 3. S1: permuted-target control

For each of the 74 targets (and the MNA diagnostic), the target vector was permuted once within the 85 lineage × medium strata and M0 + M1 were refitted under the frozen folds. Because the permutation preserves stratum means, R²(M0) stays positive (median 0.054); what must vanish is the gain from genes.

| Statistic (74 targets) | Value |
|---|---|
| median ΔR²(P) | −0.0004 |
| range | −0.0034 to +0.0084 |
| SD across targets | 0.0016 |
| median repeat-to-repeat SD within a target | 0.0008 |
| targets with ΔR²(P) > 0.02 | 0 |

The slightly positive maximum (+0.008, within noise) shows the typical magnitude of chance gains with these fold sizes, which is also why the minimum meaningful effect was set at 0.02. Full table: `results/phase2/s1_permuted_control.csv`.

## 4. S2: synthetic-target recovery

Twenty synthetic targets were built on the real gene sets of 20 targets chosen by rule (every 3rd in gene-set-size order; k from 5 to 36). Each y = lineage effect (variance 0.30) + medium effect (0.05) + a linear signal from 5 randomly chosen genes of the set (0.15) + noise (0.50). M1 used the full set, as it will in Phase 3; N-met used B = 200 (D44).

| Statistic (20 simulations) | Value |
|---|---|
| injected signal share (median) | 0.149 |
| ΔR²(M1) (median) | 0.105 |
| recovery ratio ΔR²(M1) / injected: min / median / max | 0.50 / 0.69 / 0.89 |
| correlation of recovery ratio with k | 0.00 |
| median absolute error | 0.045 (criterion ≤ 0.05) |
| simulations beating N-met at p < 0.05 | 20 / 20 |
| N-met median ΔR² under the synthetic targets | 0.010 |

**Reading.** The pipeline detects an injected local signal in every simulation and never mistakes random metabolic genes for it. It recovers about 70% of the signal rather than all of it: the ridge penalty shrinks coefficients, only 5 of the k genes carry signal, and the out-of-fold R² is estimated on 147-line test folds. This attenuation is uniform across set sizes and is the conservative behaviour the design intended (research_plan §7). It means that in Phase 3 a true ΔR² of 0.03 may be observed as ~0.02, close to the minimum meaningful effect; this will be stated as a limitation. Full table: `results/phase2/s2_synthetic.csv`.

## 5. S3: leakage checklist (research_plan §10.3)

| Item | Answer | Evidence |
|---|---|---|
| Feature selection, filtering or ranking that uses the target outside the training fold? | **No** | Gene sets come from Human-GEM; the expressed-gene filter uses expression means only; both were frozen before any outcome was read (D27, D39). |
| Standardisation, imputation or baseline fitting on all lines rather than the training fold? | **No** | `Evaluator.stats` computes gene mean/SD from training rows per fold; `target_folds` fits OLS on training rows. Test: `test_gene_standardisation_uses_training_rows_only` perturbs test rows and shows training statistics unchanged. No imputation is performed by the pipeline (the data generators' MICE imputation is a stated limitation). |
| Ridge penalty tuned on test folds? | **No** | LOO over the training fold only (`ridge.fit` receives training rows exclusively); equality with `RidgeCV(cv=None)` verified in `test_ridge_loo_matches_sklearn_ridgecv`. |
| Target list, gene-set rules, thresholds or model changed after seeing performance? | **No** | Frozen at commit `5ed608c` (D39) before any model was fitted. The only real-target result produced in Phase 2 is the S5 timing target; nothing was changed afterwards, and the deviation log is empty post-freeze. |
| Random-set nulls processed differently from M1? | **No** | `null_draws` calls the same `evaluate` with the same folds, baseline residuals and ridge grid; only the gene indices differ. |
| Duplicate or same-donor lines split across folds? | **No** | Grouped folds; `validate_folds` asserts it; duplicates were removed in Phase 0 (D16). |
| Permuted-target control shows ΔR² meaningfully above zero? | **No** | S1: median −0.0004, max 0.008. |
| Results reported for a favourable subset of targets or repeats? | **No** | Phase 3 will report all 74 targets and all repeats (expectations §7); the per-target output format stores every repeat. |

## 6. S5: timing and compute contingency

The target with the median primary gene-set size (2-deoxycytidine, k = 11) was run in full under the frozen folds (D45):

| Component | Seconds |
|---|---|
| M1 (50 fits) | 0.03 |
| M2, 12,666 genes (50 fits) | 10.4 |
| N-met, B = 1,000 (50,000 fits) | 26.1 |
| N-all, B = 1,000 (50,000 fits) | 28.3 |
| P | 0.06 |
| **Total** | **65** |

Scaling the null cost with each target's k gives a projected **1.7 h** for all 74 targets, far below the 8 h threshold, so **B stays at 1,000**. The timing was taken while the S1/S2 runs shared the CPU, so it is, if anything, an overestimate. The timing target's arrays are in `results/phase2/timing_target_2-deoxycytidine.npz`; the smoke test that preceded S5 printed its ΔR²(M1) to the console, which is disclosed here for completeness; no rule or threshold was touched afterwards.

## 7. Outputs

| File | Content |
|---|---|
| `results/phase2/folds.csv` | Frozen fold assignment (repeat, fold, CCLE_ID), 7,350 rows |
| `results/phase2/s1_permuted_control.csv`, `s1_summary.json` | Per-target permuted-target ΔR² (mean, min, max, SD over repeats); pass/fail |
| `results/phase2/s2_synthetic.csv`, `s2_summary.json` | Per-simulation injected share, ΔR²(M1), recovery ratio, p vs N-met; pass/fail |
| `results/phase2/s5_timing.json`, `timing_target_2-deoxycytidine.npz` | Timing breakdown and projection; the timing target's arrays |
| `src/phase2/{data,folds,ridge,pipeline,run}.py` | Pipeline; `run_target` is the Phase 3 entry point |
| `tests/test_phase2.py` | Ridge ≡ RidgeCV; fold partition/grouping/stratification; residual zero-mean; training-only standardisation; true genes beat random sets on toy data; within-stratum permutation |

## 8. Carried into Phase 3

1. Run `run_target` for all 74 targets and the MNA diagnostic with B = 1,000; store every array; compute ΔR², empirical p-values, BH q-values, capture ratios, the global tests (full-target primary, clustered secondary, D37), sensitivity analyses 1–4 (D36 drops 7 targets from analysis 2), and the figures (research_plan §14).
2. State the ~70% recovery (S2) as a limitation when interpreting ΔR² magnitudes.
3. The Bone lineage is absent from 3 test folds; no action needed, noted for completeness.
