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
- **Gene sets (Phase 1, D23/D24/D33):** first neighbours in Human-GEM v2.0.0 after collapsing compartments; boundary reactions and the subsystems *Exchange/demand reactions*, *Artificial reactions*, *Pool reactions* and *Protein degradation* are excluded; union of GPR genes; transporters included in the primary set. *k* counts genes that map to the expression matrix and are expressed. Frozen members: `data/processed/phase1/gene_set_members.csv`.

---

## 2. Sanity checks that must pass before any hypothesis is evaluated

| ID | Check | Pass criterion | If it fails |
|---|---|---|---|
| S1 | Permuted-target control (P) | Median ΔR² across targets within ±0.005 and no target with ΔR² > 0.02 | Leak or bug; stop and debug |
| S2 | Synthetic target recovery | A simulated target built as a noisy linear function of 5 genes in a defined set is recovered by M1 with ΔR² close to the simulated signal, and beats N-met for ≥ 95% of simulations | Pipeline cannot detect true signal; fix before proceeding |
| S3 | Leakage checklist (`research_plan.md`, Section 10.3) | All items "no" | Fix; rerun |
| S4 | Data audit | Go/no-go criteria in `research_plan.md`, Section 4 met | Narrow the question (e.g., single medium) rather than add complexity |
| S5 | Compute contingency (D38) | One complete target is timed before the rest. Projected total ≤ ~8 h at B = 1,000 | Reduce B to **500** for N-met and N-all, log the change here and in the decision log **before** running any further target |

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
- **Secondary reading (pre-specified 2026-11-16, D37):** targets are clustered by single linkage at Jaccard > 0.5 on their expressed primary gene sets (`data/processed/phase1/target_clusters.csv`: 74 targets → 56 clusters; multi-target clusters: nucleosides 9, carnitines 6, branched-chain amino acids 3, aromatic amino acids 3, asparagine/threonine 2). Each cluster contributes **one** observation: the within-cluster median of ΔR²(M1) − median ΔR²(N-met). The one-sided Wilcoxon signed-rank test is repeated on these 56 observations, and the binomial test counts clusters whose within-cluster median of ΔR²(M1) − 95th N-met percentile exceeds 0 against 5% expected. **The full-target reading stays primary** and alone determines the supported / partial / not supported call; the clustered reading is reported alongside it.

### H3 — Network genes also beat random genes from the whole transcriptome

- **Test:** As H2, against N-all.
- **Expected relation to H2:** Beating N-all is easier than beating N-met, because metabolic genes as a class carry tissue and cell-state information (Peng et al. 2018). If M1 beats N-all but not N-met, the conclusion is "metabolic genes are informative", not "network locality is informative".

### H4 — The all-gene model is at least as good as the network model

- **Rationale:** Li et al. (2022) found all-gene models generally performed better than metabolic subsets; flexible all-gene models capture lineage-correlated cell state.
- **Prediction:** ΔR²(M2) ≥ ΔR²(M1) for most targets.
- **Test:** Paired Wilcoxon signed-rank of ΔR²(M2) − ΔR²(M1) across targets; capture ratio ΔR²(M1)/ΔR²(M2) for targets with ΔR²(M2) > 0.02.
- **Expected result:** M2 better on average. The interesting outcome is the subset of targets where the capture ratio is high (e.g., ≥ 0.7), i.e., a small, interpretable gene set recovers most of the achievable signal.
- **Surprising (and to be scrutinised for leakage) if:** M1 substantially beats M2 for many targets.

### H5 — Literature positive controls: does the first-neighbour set capture a known association?

*Amended 2026-11-16 (D34), before any modelling; the original wording is in the deviation log. In Human-GEM, three of the five literature genes are not first neighbours of their metabolite, so H5 now names the gene(s) that are in each set.*

- **Controls (all pass selection rules 1–4):**
  - **kynurenine** — in-set genes **AFMID, KYNU** (IDO1 and TDO2 are at network distance 2: they produce N-formylkynurenine, which AFMID hydrolyses to kynurenine);
  - **creatine** — in-set genes **GAMT, SLC6A8** (GATM is at distance 2: it makes guanidinoacetate);
  - **taurine** — in-set gene **SLC6A6**.
- **Excluded by rule 3 (not controls):** 1-methylnicotinamide (gene set NNMT, AOX1: 2 expressed genes) and N-carbamoyl-β-alanine (DPYS and UPB1 in set but not expressed; DPYD at distance 2).
- **Prediction:** each of the three controls has ΔR²(M1) ≥ 0.02 and empirical p vs N-met < 0.05 (uncorrected; confirmatory).
- **If controls fail:** either the pipeline is too conservative (check S2) or the two-stage lineage/medium adjustment absorbs the signal. Report both M0-adjusted and unadjusted results for controls to diagnose which.
- **Pipeline diagnostic (not a target, outside selection rules):** M1 is fitted for **1-methylnicotinamide** with its 2-gene set (NNMT, AOX1) under the identical pipeline, and ΔR²(M1) with its N-met null is reported. It is labelled as a diagnostic, enters none of H1–H4 and is not counted in any target tally.
- **Stretch goal 2 motivation:** kynurenine's literature genes being at distance 2 is the concrete case for second-neighbour gene sets (`research_plan.md` §12).

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
2. Enzymatic-only gene sets (transport reactions excluded). **Targets with < 3 expressed enzymatic genes are dropped from this analysis only (D36):** carnosine, glucuronate, histidine, oxalate, taurine, thymine, tryptophan. 67 of 74 targets are analysed.
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
| 2026-11-16 | Gene-set rule: *Protein degradation* subsystem added to the exclusions (D33) | 35 lumped plasma-protein hydrolysis reactions attached the same proteases to all amino-acid targets | **Before** (no model fitted) |
| 2026-11-16 | H5 reworded (D34): controls are tested on their in-set genes (kynurenine → AFMID/KYNU; creatine → GAMT/SLC6A8; taurine → SLC6A6); MNA and N-carbamoyl-β-alanine excluded by rule 3; MNA kept as a labelled pipeline diagnostic outside H1–H4. Original: "Targets (if they pass selection rules): 1-methylnicotinamide (NNMT), kynurenine (IDO1/TDO2), creatine (SLC6A8/GATM), N-carbamoyl-β-alanine (DPYD), taurine (SLC6A6)." | Three literature genes are not first neighbours in Human-GEM | **Before** (no model fitted) |
| 2026-11-16 | H2: pre-specified secondary clustered reading added (D37); full-target reading stays primary | Overlapping gene sets | **Before** (no model fitted) |
| 2026-11-16 | Sensitivity analysis 2: targets with < 3 expressed enzymatic genes dropped from that analysis (D36) | Degenerate size-matched null | **Before** (no model fitted) |
| 2026-11-16 | Compute contingency S5 added (D38): B may be reduced to 500 after timing one target | Keep Phase 2 finishable; decided before any run | **Before** (no model fitted) |

---

## 9. Target list (to be appended at freeze)

**Frozen on: 2026-11-16 (commit: `5ed608cd126467a3df41810c37204c2f86ebeaf7`).** Generated from `data/processed/phase1/target_selection.csv` by `src/phase1/report_tables.py` after the project owner's review decisions D31–D38. **74 targets.** From this point every change to the list is a logged deviation (§8).

- *k* = expressed genes in the first-neighbour set (primary / enzymatic-only); "S2" marks targets kept in sensitivity analysis 2 (D36).
- Cluster = single-linkage Jaccard > 0.5 cluster for the secondary H2 reading (D37); blank = singleton.
- Positive controls (H5): kynurenine, creatine, taurine. Pipeline diagnostic, not a target: 1-methylnicotinamide (NNMT, AOX1).

<!-- TARGETS:START -->
| # | Target (CCLE feature) | HMDB | KEGG | Human-GEM ID | Human-GEM name | k (primary) | k (enzymatic) | S2 | Class | Pos. control | Cluster |
|---:|---|---|---|---|---|---:|---:|---|---|---|---|
| 1 | choline | HMDB0000097 | C00114 | MAM01513 | choline | 15 | 9 | yes | medium-supplied |  |  |
| 2 | histidine | HMDB0000177 | C00135 | MAM02125 | histidine | 14 | 2 | dropped | medium-supplied |  |  |
| 3 | isoleucine | HMDB0000172 | C00407 | MAM02184 | isoleucine | 12 | 3 | yes | medium-supplied |  | C03 |
| 4 | leucine | HMDB0000687 | C00123 | MAM02360 | leucine | 16 | 5 | yes | medium-supplied |  | C03 |
| 5 | lysine | HMDB0000182 | C00047 | MAM02426 | lysine | 16 | 7 | yes | medium-supplied |  |  |
| 6 | methionine | HMDB0000696 | C00073 | MAM02471 | methionine | 26 | 14 | yes | medium-supplied |  |  |
| 7 | pantothenate | HMDB0000210 | C00864 | MAM02680 | pantothenate | 4 | 3 | yes | medium-supplied |  |  |
| 8 | phenylalanine | HMDB0000159 | C00079 | MAM02724 | phenylalanine | 15 | 5 | yes | medium-supplied |  | C04 |
| 9 | threonine | HMDB0000167 | C00188 | MAM02993 | threonine | 15 | 3 | yes | medium-supplied |  | C05 |
| 10 | tryptophan | HMDB0000929 | C00078 | MAM03089 | tryptophan | 12 | 2 | dropped | medium-supplied |  | C04 |
| 11 | valine | HMDB0000883 | C00183 | MAM03135 | valine | 12 | 3 | yes | medium-supplied |  | C03 |
| 12 | arginine | HMDB0000517 | C00062 | MAM01365 | arginine | 13 | 6 | yes | mixed |  |  |
| 13 | asparagine | HMDB0000168 | C00152 | MAM01369 | asparagine | 17 | 4 | yes | mixed |  | C05 |
| 14 | aspartate | HMDB0000191 | C00049 | MAM01370 | aspartate | 17 | 12 | yes | mixed |  |  |
| 15 | glutamate | HMDB0000148 | C00025 | MAM01974 | glutamate | 58 | 50 | yes | mixed |  |  |
| 16 | glutamine | HMDB0000641 | C00064 | MAM01975 | glutamine | 30 | 17 | yes | mixed |  |  |
| 17 | glutathione reduced | HMDB0000125 | C00051 | MAM02026 | GSH | 53 | 49 | yes | mixed |  |  |
| 18 | glycine | HMDB0000123 | C00037 | MAM01986 | glycine | 38 | 27 | yes | mixed |  |  |
| 19 | inositol | HMDB0000211 | C00137 | MAM02171 | inositol | 8 | 6 | yes | mixed |  |  |
| 20 | niacinamide | HMDB0001406 | C00153 | MAM02583 | nicotinamide | 24 | 24 | yes | mixed |  |  |
| 21 | proline | HMDB0000162 | C00148 | MAM02770 | proline | 20 | 9 | yes | mixed |  |  |
| 22 | serine | HMDB0000187 | C00065 | MAM02896 | serine | 29 | 14 | yes | mixed |  |  |
| 23 | tyrosine | HMDB0000158 | C00082 | MAM03101 | tyrosine | 15 | 7 | yes | mixed |  | C04 |
| 24 | 2-deoxyadenosine | HMDB0000101 | C00559 | MAM01666 | deoxyadenosine | 10 | 8 | yes | synthesised |  | C01 |
| 25 | 2-deoxycytidine | HMDB0000014 | C00881 | MAM01668 | deoxycytidine | 11 | 9 | yes | synthesised |  | C01 |
| 26 | 3-phosphoglycerate | — | C00197 | MAM00913 | 3-phospho-D-glycerate | 9 | 9 | yes | synthesised |  |  |
| 27 | 5-adenosylhomocysteine | HMDB0000939 | C00021 | MAM02871 | SAH | 59 | 58 | yes | synthesised |  |  |
| 28 | CMP | HMDB0000095 | C00055 | MAM01590 | CMP | 36 | 34 | yes | synthesised |  |  |
| 29 | GABA | HMDB0000112 | C00334 | MAM00970 | 4-aminobutyrate | 9 | 8 | yes | synthesised |  |  |
| 30 | GMP | HMDB0001397 | C00144 | MAM02016 | GMP | 34 | 32 | yes | synthesised |  |  |
| 31 | PEP | HMDB0000263 | C00074 | MAM02696 | PEP | 7 | 6 | yes | synthesised |  |  |
| 32 | acetylcarnitine | HMDB0000201 | C02571 | MAM02634 | O-acetylcarnitine | 7 | 5 | yes | synthesised |  | C02 |
| 33 | adenosine | HMDB0000050 | C00212 | MAM01280 | adenosine | 18 | 14 | yes | synthesised |  | C01 |
| 34 | alanine | HMDB0000161 | C00041 | MAM01307 | alanine | 36 | 22 | yes | synthesised |  |  |
| 35 | alpha-glycerophosphate | — | C00093 | MAM02914 | sn-glycerol-3-phosphate | 11 | 11 | yes | synthesised |  |  |
| 36 | alpha-glycerophosphocholine | HMDB0000086 | C00670 | MAM02912 | sn-glycerol-3-PC | 8 | 8 | yes | synthesised |  |  |
| 37 | alpha-ketoglutarate | HMDB0000208 | C00026 | MAM01306 | AKG | 47 | 36 | yes | synthesised |  |  |
| 38 | carnitine | HMDB0000062 | C00318 | MAM02348 | L-carnitine | 13 | 6 | yes | synthesised |  | C02 |
| 39 | carnosine | HMDB0000033 | C00386 | MAM01423 | carnosine | 4 | 2 | dropped | synthesised |  |  |
| 40 | citrate | HMDB0000094 | C00158 | MAM01587 | citrate | 6 | 4 | yes | synthesised |  |  |
| 41 | citrulline | HMDB0000904 | C00327 | MAM01588 | citrulline | 7 | 3 | yes | synthesised |  |  |
| 42 | creatine | HMDB0000064 | C00300 | MAM01619 | creatine | 6 | 4 | yes | synthesised | yes |  |
| 43 | cystathionine | HMDB0000099 | C02291 | MAM02349 | L-cystathionine | 3 | 3 | yes | synthesised |  |  |
| 44 | cytidine | HMDB0000089 | C00475 | MAM01630 | cytidine | 12 | 9 | yes | synthesised |  | C01 |
| 45 | erythrose-4-phosphate | — | C00279 | MAM01785 | erythrose-4-phosphate | 3 | 3 | yes | synthesised |  |  |
| 46 | glucuronate | HMDB0000127 | C00191 | MAM01973 | glucuronate | 3 | 2 | dropped | synthesised |  |  |
| 47 | glutathione oxidized | HMDB0003337 | C00127 | MAM02027 | GSSG | 18 | 18 | yes | synthesised |  |  |
| 48 | guanosine | HMDB0000133 | C00387 | MAM02038 | guanosine | 11 | 8 | yes | synthesised |  | C01 |
| 49 | hexanoylcarnitine | HMDB0000705 (curator, unverified) | — | MAM03498 | Hexanoyl Carnitine | 5 | 4 | yes | synthesised |  |  |
| 50 | homocysteine | HMDB0000742 | C00155 | MAM02133 | homocysteine | 6 | 6 | yes | synthesised |  |  |
| 51 | hypoxanthine | HMDB0000157 | C00262 | MAM02159 | hypoxanthine | 5 | 3 | yes | synthesised |  |  |
| 52 | inosine | HMDB0000195 | C00294 | MAM02170 | inosine | 9 | 6 | yes | synthesised |  | C01 |
| 53 | isocitrate | HMDB0000193 | C00311 | MAM02183 | isocitrate | 9 | 7 | yes | synthesised |  |  |
| 54 | kynurenic acid | HMDB0000715 | C01717 | MAM00990 | 4-hydroxy-2-quinolinecarboxylic acid | 3 | 3 | yes | synthesised |  |  |
| 55 | kynurenine | HMDB0000684 | C00328 | MAM02319 | kynurenine | 7 | 5 | yes | synthesised | yes |  |
| 56 | lactate | HMDB0000190 | C00186 | MAM02403 | L-lactate | 16 | 10 | yes | synthesised |  |  |
| 57 | malate | HMDB0000156 | C00149 | MAM02439 | malate | 11 | 8 | yes | synthesised |  |  |
| 58 | malonylcarnitine | HMDB0002095 (curator, unverified) | — | MAM02443 | malonyl-carnitine | 7 | 5 | yes | synthesised |  | C02 |
| 59 | ornithine | HMDB0000214 | C00077 | MAM02658 | ornithine | 10 | 5 | yes | synthesised |  |  |
| 60 | oxalate | HMDB0002329 | C00209 | MAM02661 | oxalate | 6 | 1 | dropped | synthesised |  |  |
| 61 | palmitoylcarnitine | HMDB0000222 | C02990 | MAM02411 | L-palmitoylcarnitine | 8 | 5 | yes | synthesised |  | C02 |
| 62 | phosphocreatine | HMDB0001511 | C02305 | MAM01620 | creatine-phosphate | 3 | 3 | yes | synthesised |  |  |
| 63 | propionylcarnitine | HMDB0000824 | C03017 | MAM02657 | O-propanoylcarnitine | 7 | 5 | yes | synthesised |  | C02 |
| 64 | putrescine | HMDB0001414 | C00134 | MAM02812 | putrescine | 10 | 8 | yes | synthesised |  |  |
| 65 | sorbitol | HMDB0000247 | C00794 | MAM01682 | D-glucitol | 7 | 3 | yes | synthesised |  |  |
| 66 | stearoylcarnitine | HMDB0000848 | — | MAM02940 | stearoylcarnitine | 8 | 5 | yes | synthesised |  | C02 |
| 67 | taurine | HMDB0000251 | C00245 | MAM02961 | taurine | 3 | 1 | dropped | synthesised | yes |  |
| 68 | thymidine | HMDB0000273 | C00214 | MAM02996 | thymidine | 13 | 10 | yes | synthesised |  | C01 |
| 69 | thymine | HMDB0000262 (HG) | C00178 | MAM02997 | thymine | 3 | 2 | dropped | synthesised |  |  |
| 70 | uracil | HMDB0000300 | C00106 | MAM03118 | uracil | 5 | 4 | yes | synthesised |  |  |
| 71 | urate | HMDB0000289 | C00366 | MAM03120 | urate | 5 | 3 | yes | synthesised |  |  |
| 72 | uridine | HMDB0000296 | C00299 | MAM03123 | uridine | 14 | 11 | yes | synthesised |  | C01 |
| 73 | xanthine | HMDB0000292 | C00385 | MAM03148 | xanthine | 4 | 4 | yes | synthesised |  |  |
| 74 | xanthosine | HMDB0000299 | C01762 | MAM03149 | xanthosine | 6 | 6 | yes | synthesised |  | C01 |
<!-- TARGETS:END -->
