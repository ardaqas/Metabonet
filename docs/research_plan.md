# Research Plan

**Project:** Predicting metabolite abundance from gene expression using human metabolic network information (CCLE)
**Status:** Draft v0.1 (16 November 2026). Nothing has been implemented, downloaded or run.
**Companion documents:** `literature_review.md` (evidence base), `expectations.md` (pre-specified hypotheses and decision criteria)

---

## 0. Summary and central recommendation

**The smallest scientifically meaningful version of this project:**

> For roughly 20–35 **polar, unambiguously identified CCLE metabolites** that map one-to-one to a pinned Human-GEM release, test whether the expression of genes in reactions that produce, consume or transport each metabolite improves **held-out** prediction beyond **lineage and culture medium**, and whether that improvement exceeds what **size-matched random gene sets** (random metabolic genes and random expressed genes) achieve under identical cross-validation. Compare against an **all-gene** model as a ceiling. Use **one** model class (ridge regression), **one** dataset pair (CCLE 2019 metabolomics + matched RNA-seq), and **solid-tumour cell lines**, with a permuted-target negative control and a handful of literature positive controls.

Why this is the right size:
- It answers one clear question that prior work has not directly answered (`literature_review.md`, Section 9.3), using public data and laptop-scale compute.
- Every component demonstrates a skill reviewers look for: biological curation (metabolite mapping, GPRs, target selection), data engineering (identifier matching), statistics (nested CV, empirical nulls, multiple testing), and critical interpretation (likely mixed or negative results).
- Nothing in it depends on a fragile or ambitious step (no flux simulation, no context-specific models, no deep learning, no external validation cohort).

What to exclude, and why, is in Section 12.

---

## 1. Research question

**Primary question.** For a pre-specified set of metabolites, does the expression of each metabolite's first-neighbour network genes (from Human-GEM) carry predictive information about its abundance that is (a) additional to lineage and culture medium and (b) greater than the information carried by equally sized random gene sets?

**Secondary questions.**
1. How does the network gene set compare with an all-gene model (ceiling), i.e., what share of the achievable expression-based improvement is captured by a small, interpretable gene set?
2. Which metabolite properties are associated with network-specific predictability (e.g., dedicated enzyme vs hub metabolite; medium-supplied nutrient vs synthesised intermediate)?

**Not asked.** Whether expression causes metabolite changes; whether results transfer to tumours; whether flux can be inferred.

---

## 2. Feasibility of combining CCLE metabolomics and transcriptomics

### 2.1 Verdict
**Feasible for the proposed scope, with known and manageable limitations.** The data are public, matched on the same cell line collection, large enough (several hundred solid-tumour lines) for per-metabolite regression with ≤ 60 features, and have already been combined in at least three published analyses.

### 2.2 Strengths
- Same CCLE cell line collection and same project era for both data types (metabolomics: Li et al. 2019; RNA-seq: Ghandi et al. 2019, 1,019 lines).
- Targeted LC-MS panel with standards-based identification and good technical reproducibility (221/225 metabolites with CV < 0.1 in replicated lines).
- Polar metabolites (124) cover central carbon, amino acid, nucleotide and cofactor metabolism, which are well represented in Human-GEM.

### 2.3 Limitations that shape the design

| Limitation | Consequence | Design response |
|---|---|---|
| Metabolomics and RNA-seq come from different cultures and time points | Adds noise; caps achievable R² | Accept; interpret absolute R² conservatively; focus on *relative* comparisons between models |
| Relative (normalised) abundances | Cannot compare across metabolites or convert to concentrations | Model each metabolite separately; compare models within a metabolite |
| Lines grown in ~18 different media | Medium confounds expression–metabolite links | Medium category as a baseline covariate; medium-supplied nutrients interpreted separately |
| Medium is collinear with lineage | Effects hard to separate | Model them jointly in the baseline; do not interpret their individual coefficients |
| ~9% values imputed (MICE) by data generators | Small information sharing across samples; imputed values less reliable | Cannot undo; exclude targets with high missingness if identifiable from supplementary data |
| Four LC-MS batches; batch labels may not be public | Possible residual batch structure | Check whether batch is available in supplementary data; if so, include as covariate in a sensitivity analysis. **Phase 0 audit: not available; analysis dropped and stated as a limitation (D19)** |
| Suspension (mostly haematopoietic) vs adherent growth | Systematic metabolic and culture differences (Li et al. 2022 excluded these) | Primary analysis on solid-tumour lines; all-lines analysis as sensitivity |
| Polar/lipid counts differ between papers (124/101 vs 136/89) | Risk of misclassifying targets | Classify by assay method in the original supplementary table |
| Possible duplicate, misidentified or same-donor lines | Leakage across folds; inflated performance | Check Cellosaurus problematic-line flags and shared-donor annotations; keep same-donor lines in the same fold |

### 2.4 Recommended data sources (to be confirmed in Phase 0; **not downloaded yet**)
- **Metabolomics:** DepMap portal, "CCLE 2019" release, `CCLE_metabolomics_20190502.csv`; plus Li et al. (2019) Supplementary Table 1 (assay method per metabolite, medium per line, CVs).
- **Expression:** CCLE 2019 RNA-seq, RSEM gene-level TPM (1,019 lines; Ensembl gene IDs, GENCODE v19). *Rationale:* same release and identifier namespace as the metabolomics. *Trade-off:* older annotation than current DepMap releases; mapping to current Ensembl IDs used by Human-GEM will lose a small number of genes. Using a current DepMap expression release is a defensible alternative; whichever is chosen, record it in the decision log and do not mix.
- **Cell line annotation:** `Cell_lines_annotations_20181226.txt` (CCLE_ID ↔ depMapID, `Growth.Medium`, `Supplements`, site/histology) and a current DepMap `Model.csv` (ModelID, CCLEName, OncotreeLineage, RRID).
- **Network:** Human-GEM, one tagged release (recommend **v2.0.0**, the most curated; record the commit hash), SBML file plus `metabolites.tsv` and `genes.tsv`.

---

## 3. Are Human-GEM and COBRApy appropriate?

### 3.1 Verdict
**Yes, with explicit rules.** Human-GEM is the most carefully curated open human GEM, it has external metabolite identifiers needed for mapping, it uses Ensembl gene IDs, and it was used for exactly this purpose (gene–metabolite network distance) by Benedetti et al. (2023). COBRApy is a standard, well-documented way to load the SBML file and query metabolite → reaction → gene relationships. No optimisation or solver is needed.

### 3.2 What COBRApy will be used for (and not used for)
- **Used for:** loading the SBML model; iterating over each metabolite's reactions; reading each reaction's gene set and GPR rule; identifying exchange, sink, demand and biomass/pseudo-reactions to exclude.
- **Not used for:** FBA, FVA, gene deletions, context-specific model building, or any simulation.
- Pandas will read `metabolites.tsv` (external IDs) and `genes.tsv` (gene symbols) alongside.

### 3.3 Known problems and the rule adopted for each

| Issue | Why it matters | Rule |
|---|---|---|
| **Compartment-specific metabolite copies** (e.g., the same metabolite in cytosol, mitochondria, extracellular space) | CCLE measures whole-cell extracts, not compartments | Collapse all compartment copies of a metabolite to one node; take the union of their reactions |
| **Transport reactions** | After collapsing, transport looks like X → X; but transporters are biologically relevant (e.g., SLC6A8–creatine, SLC6A6–taurine) | **Primary:** include transporter genes. **Sensitivity:** enzymatic reactions only |
| **Exchange / sink / demand / biomass / pseudo-reactions** | Modelling artefacts, not biology | Exclude. *Phase 1 (D23, D33): implemented as boundary reactions plus the subsystems Exchange/demand, Artificial, Pool and **Protein degradation** (lumped plasma-protein hydrolysis)* |
| **GPR logic (AND = complex subunits; OR = isozymes)** | Feature-set construction ignores logic; large complexes inflate sets | Take the union of all genes in the GPR. Report set sizes. Reaction-level scores (min over AND, max over OR) are **out of scope** for the primary analysis |
| **Reactions without GPR** (spontaneous or orphan) | Contribute no genes | Allowed; record how many reactions per target lack genes |
| **Promiscuous genes** (generic transporters, broad-specificity enzymes) | Same genes appear in many targets' sets; targets are not independent | Report pairwise overlap (Jaccard) between target sets; interpret the global test with this in mind |
| **Hub / currency metabolites** (ATP, NAD(P)(H), CoA, H₂O, Pi, CO₂, H⁺, O₂; also glutamate, 2-oxoglutarate in transaminations) | Gene sets become enormous and non-specific | Exclude as **targets** via a pre-specified gene-set size cap (Section 5). Because only first neighbours are used, currency metabolites do not need to be removed from the graph |
| **GPR errors and version drift** (Luo et al. 2026: 1,135 GPRs updated in v2.0.0) | Gene sets depend on version | Pin one release; optional check of set stability between v1.19 and v2.0 |
| **Gene identifier mapping** (Human-GEM current Ensembl IDs vs GENCODE v19 IDs in CCLE 2019) | Silent loss of genes | Strip version suffixes; map via Ensembl ID, then fall back to HGNC symbol; report mapping rate per target and overall |
| **Direction (production vs consumption)** | Many reactions are reversible; direction is not identifiable from expression | Do not split by direction |

### 3.4 Metabolite mapping (the highest-risk step)
CCLE metabolite names must be mapped to Human-GEM metabolites by hand-checked identifiers, not string matching.

Procedure (manual, documented in a mapping table):
1. For each polar CCLE metabolite, assign HMDB, KEGG and ChEBI identifiers from the original supplementary data where available, otherwise by manual lookup.
2. Match to Human-GEM through `metabolites.tsv` external identifiers (HMDB first, then KEGG/ChEBI), at the compartment-free metabolite level.
3. Manually verify each match: stereochemistry (L- vs D-), generic vs specific entries (e.g., "hexose" vs glucose), protonation state (irrelevant), and whether the LC-MS feature could represent co-eluting isomers.
4. Assign each metabolite a mapping status: **1:1 confident**, **ambiguous** (isomers or one-to-many), or **unmapped**. Only 1:1 confident mappings become targets.
5. Record the rationale for every exclusion. This table is itself a portfolio deliverable.

---

## 4. Data matching protocol (Phase 0 audit)

Goal: a single, reproducible table of cell lines with metabolomics, expression, lineage and medium, before any modelling.

1. Load metabolomics; record number of lines, metabolites, missing values, and value distribution (check whether values are already log-transformed; do **not** log twice).
2. Map metabolomics `CCLE_ID` → DepMap ID using the 2019 annotation file; cross-check with current `Model.csv` (`CCLEName`). Log all unmatched or conflicting IDs.
3. Intersect with expression; record the overlap.
4. Assign lineage (OncotreeLineage via DepMap ID). Collapse lineages with fewer than 15 lines into "other".
5. Assign medium: prefer Li et al. (2019) Supplementary Table 1; fall back to `Growth.Medium`. Normalise free text to categories (e.g., RPMI-1640, DMEM, EMEM/MEM, McCoy's 5A, DMEM:F12, Leibovitz L-15, IMDM, Williams' E, other). Collapse categories with fewer than 15 lines into "other". Report agreement between the two sources.
6. Flag growth type (adherent/suspension) and haematopoietic/lymphoid lineage.
7. Flag duplicates, problematic lines (Cellosaurus) and shared-donor groups.
8. Freeze the cell line table (version and hash).

**Go/no-go checkpoint.** Proceed only if: ≥ 400 solid-tumour lines have metabolomics, expression, lineage and medium; medium can be assigned for ≥ 90% of them; and ≥ 15 polar metabolites pass target selection (Section 5). If not, narrow the question (e.g., single largest medium, RPMI only) rather than adding complexity.

---

## 5. Target selection (pre-specified; outcome-blind)

A polar metabolite becomes a target only if **all** apply:
1. Measured by a HILIC (polar) method.
2. Mapping status **1:1 confident** (Section 3.4); not an isomer-ambiguous or summed feature.
3. After collapsing compartments and applying the rules in Section 3.3, its first-neighbour gene set contains **3 to 60 genes** present in the expression matrix. (Fewer than 3 makes the random-set null degenerate; more than 60 indicates a hub metabolite with non-specific neighbourhood.)
4. Missingness acceptable and non-trivial variance across lines. **Specified 2026-11-16 (D18): ≤ 20 % of primary-set values missing before MICE imputation** (`expectations.md` §1a).

These rules use no information about predictive performance. Expected yield: roughly 20–35 targets; the exact number is whatever the rules produce. **No targets may be added or removed after modelling starts**; any later change is reported as a deviation.

Each target is annotated (before modelling) with:
- **Class:** synthesised intermediate; medium-supplied nutrient (e.g., essential amino acids); or mixed.
- **Literature positive control:** yes/no (MNA–NNMT; kynurenine–IDO1/TDO2; creatine–SLC6A8/GATM; N-carbamoyl-β-alanine–DPYD; taurine–SLC6A6), conditional on passing the rules above.

---

## 6. Feature sets

For each target metabolite *m* with network gene set *Gₘ* of size *k*:

| Label | Features | Purpose |
|---|---|---|
| **M0** baseline | Lineage + medium (one-hot) | What is explained without expression |
| **M1** network | M0 + expression of *Gₘ* | Hypothesis model |
| **M2** all-gene | M0 + expression of all expressed protein-coding genes | Ceiling for expression-based prediction |
| **N-met** null | M0 + *k* random genes drawn from Human-GEM metabolic genes (excluding *Gₘ*) | Tests network locality vs "metabolic genes in general" (the key control; see Peng 2018, Recon8D) |
| **N-all** null | M0 + *k* random genes from all expressed genes (excluding *Gₘ*) | Tests whether *Gₘ* beats arbitrary genes |
| **P** negative control | M1 fitted to the target permuted within lineage × medium strata | Checks the pipeline cannot manufacture signal |

Expression: log2(TPM + 1), protein-coding genes, with the 47 GENCODE v19 `ENSGR` chrY pseudoautosomal rows dropped first (D22); "expressed" = mean log2(TPM+1) ≥ 1 across lines (an outcome-blind filter; fixed before modelling). Random sets: B = 1,000 per target for N-met and N-all (compute is trivial at this size). Sensitivity: variance-matched sampling (draw random genes from the same expression-variance deciles as *Gₘ*).

---

## 7. Model

**One model class: ridge regression**, with a two-stage structure that makes "beyond lineage and medium" explicit:

Within each training fold:
1. Fit the baseline (ordinary least squares on one-hot lineage + medium) to the target.
2. Compute training residuals.
3. Standardise gene features using training-fold statistics only.
4. Fit ridge regression of the residuals on the gene features; choose the penalty by inner cross-validation (e.g., `RidgeCV` on the training fold only).
5. Predict the test fold as: baseline prediction + ridge prediction.

Why ridge: stable with correlated predictors (complex subunits, isozymes), uses all genes rather than selecting a few (avoids the LASSO instability seen in Li 2022), fast enough to run thousands of null models, and easy to explain. Why two-stage: genes cannot "steal" variance that lineage or medium already explain, so the comparison is conservative and interpretable. The same procedure is used for M1, M2, N-met, N-all and P.

**Explicitly not used:** random forests, gradient boosting, neural networks, feature selection by univariable correlation.

---

## 8. Validation design

- **Repeated K-fold cross-validation:** 5 folds × 10 repeats, stratified by lineage. Same-donor or duplicate lines are kept within one fold (grouped splitting).
- **Identical folds for every model and every random set** within a repeat, so comparisons are paired.
- **Nested tuning:** ridge penalty chosen inside each training fold only.
- **Everything learned from data happens inside the training fold:** baseline coefficients, standardisation, penalty choice. The only pre-fold steps are outcome-blind (expression filter, category collapsing, target selection), and they are frozen before modelling.
- **No separate held-out test set** is needed given repeated CV and the absence of model selection across model families; the target list, feature definitions and model are fixed in advance (`expectations.md`).

Optional (stretch, Section 12): leave-one-lineage-out, to test whether network-gene information transfers to unseen lineages (the baseline then cannot use lineage).

---

## 9. Metrics and statistical tests

**Per target, per model:** out-of-fold R² (pooled across folds within a repeat, then averaged across repeats) and Spearman ρ as a secondary metric.

**Primary quantity:** incremental R² over baseline, **ΔR²(X) = R²(M0 + X) − R²(M0)**.

**Per-target tests:**
- *Network vs null:* empirical p-value = (1 + number of random sets with ΔR² ≥ observed ΔR²(M1)) / (B + 1), separately against N-met and N-all. Benjamini–Hochberg FDR across targets.
- *Network vs ceiling:* "capture ratio" = ΔR²(M1) / ΔR²(M2), reported only where ΔR²(M2) > 0.02.

**Global tests (across targets):**
- One-sided Wilcoxon signed-rank test of ΔR²(M1) − median ΔR²(N-met) across targets.
- Number of targets exceeding the 95th percentile of their own N-met null, compared with the 5% expected under the null (binomial test). Note that overlapping gene sets make targets not fully independent; report overlap alongside.

**Stability:** spread of ΔR² across the 10 repeats (reported as an interval of repeat-to-repeat variation; not interpreted as a formal confidence interval).

Thresholds and decision rules are fixed in `expectations.md`.

---

## 10. Major risks

### 10.1 Scientific risks

| Risk | Likelihood | Consequence | Mitigation |
|---|---|---|---|
| Most predictive signal is general cell state, not local network | High (Recon8D random-gene result; Benedetti 2023) | M1 ≈ N-met for most targets | This *is* a valid answer; N-met null makes it detectable rather than misread |
| Lineage/medium explain most variance; little left for genes | Medium–high | Small ΔR² values everywhere | Report M0 R² per target; interpret ΔR² relative to residual variance |
| Steady-state pool size is decoupled from enzyme level (distributed flux control, allostery, enzyme excess) | High in principle | Weak relationships even for "correct" genes | Expected; positive controls show where single-enzyme control exists |
| Medium-supplied nutrients reflect supply and uptake | High for essential amino acids | Network genes uninformative | Pre-annotated class; analysed separately |
| Genetic/epigenetic drivers not visible in expression | Medium | Misses (e.g., 2HG–IDH) | Interpret, do not "fix" by adding mutation data in the MVP |
| mRNA–protein discordance | Medium | Attenuates signal | Out of scope; mention as limitation |

### 10.2 Technical risks

| Risk | Consequence | Mitigation / detection |
|---|---|---|
| **Cell line mismatch** (name formats, aliases, re-annotated lines) | Wrong pairings; silent noise | ID crosswalk via DepMap ID; log every unmatched/conflicting name; spot-check 20 lines manually |
| **Metabolite misidentification or isomers** | Wrong gene set for the measured molecule | 1:1 confident mapping only; isomer-ambiguous features excluded; documented mapping table |
| **Gene ID loss** (version drift) | Incomplete gene sets | Report mapping rate; flag targets losing > 20% of genes |
| **GPR artefacts** (complex inflation, promiscuous genes, gap-filled reactions) | Non-specific sets | Size cap 3–60; overlap report; enzymatic-only sensitivity; pinned version |
| **Double log-transformation or scale errors** | Distorted targets | Inspect distributions in Phase 0 |
| **Unbalanced or tiny lineage/medium categories** | Unstable baseline | Collapse categories < 15 lines; stratified folds |
| **Scope creep** | Unfinished project | Section 12 exclusions; stretch goals only after all deliverables exist |

### 10.3 Data leakage checklist (must all be "no" before results are trusted)

- [ ] Any feature selection, filtering or ranking that uses the target outside the training fold? (Li 2022's univariable pre-filter is the cautionary example.)
- [ ] Standardisation, imputation or baseline fitting done on all lines rather than the training fold?
- [ ] Ridge penalty tuned on test folds?
- [ ] Target list, gene-set rules, thresholds or model changed after seeing performance? (Any change must be logged as a deviation and reported.)
- [ ] Random-set nulls processed differently from M1 (different folds, tuning or scaling)?
- [ ] Duplicate or same-donor lines split across training and test folds?
- [ ] Permuted-target control (P) shows ΔR² meaningfully above zero? (If yes, there is a leak; stop and debug.)
- [ ] Results reported only for a favourable subset of targets or repeats?

---

## 11. Biological interpretation plan

After the statistics, select **2–3 case studies** (chosen by pre-specified rule: the top target, a positive control, and one clear failure), and for each:
- Show the gene set, which genes carry the largest ridge coefficients, and whether their direction makes biochemical sense (e.g., a producing enzyme positively and a consuming enzyme negatively associated).
- Relate to known biology (e.g., NNMT methylates nicotinamide using SAM; IDO1/TDO2 initiate tryptophan → kynurenine).
- For a failure, give the most likely reason (medium supply, genetic driver, isomer ambiguity, hub metabolite, post-transcriptional control).

Interpretation must stay associative. Coefficient signs in correlated ridge models are suggestive, not causal.

---

## 12. Scope: included, excluded, stretch

### Included (minimum viable project)
- CCLE 2019 metabolomics (polar panel only) + matched RNA-seq; solid-tumour lines.
- Human-GEM (one pinned release) parsed with COBRApy; first-neighbour gene sets.
- Ridge regression; M0, M1, M2; N-met and N-all random-set nulls; permuted-target control.
- Repeated, stratified, nested CV with identical folds.
- Pre-specified hypotheses and thresholds (`expectations.md`).
- Sensitivity analyses limited to four: (1) all lines including haematopoietic; (2) enzymatic-only gene sets; (3) variance-matched random sets; (4) excluding Cellosaurus-problematic lines (added 2026-11-16 before any results, D20).
- 2–3 biological case studies.

### Excluded (and why)

| Excluded | Reason |
|---|---|
| Lipids (101 species) | Mapping to Human-GEM lipid pseudo-species is ambiguous; lipid GPRs were heavily revised in v2.0 (717 changes); doubles curation effort |
| Flux balance analysis, context-specific GEMs (GIMME, iMAT, INIT), flux sampling | Different question (flux, not pool size); many arbitrary choices; would dominate the project |
| Random forests, boosting, deep learning, graph neural networks | Do not address the question; harder to control and explain; Recon8D already explored flexible models |
| Other omics (proteomics, methylation, CNV, mutations) | Scope; changes the question from "network genes" to "multi-omics" |
| External datasets (NCI-60, Cherkaoui 2022) | Different platforms and panels; low overlap; FIA-MS identification at confidence level 4 |
| TCGA or tissue application | No measured metabolites to validate against |
| Second-neighbour or longer network distances | Requires currency-metabolite handling and multiplies choices; stretch only |
| Reaction-level GPR scoring (min/max rules) | Adds modelling choices; gene-level features suffice for the primary question |
| Directional (production vs consumption) analysis | Reversibility makes it unidentifiable from expression |

### Stretch goals (only after all core deliverables exist; pick at most one)
1. Leave-one-lineage-out evaluation.
2. Second-neighbour gene sets (distance 2, after removing currency metabolites by a stated list). *Motivation from Phase 1 (D34): kynurenine's literature genes IDO1/TDO2 are at distance 2 in Human-GEM.*
3. Gene-set stability across Human-GEM v1.19 vs v2.0 and its effect on ΔR².

---

## 13. Work plan

Assumes ~10 hours per week.

| Phase | Weeks | Output | Checkpoint |
|---|---|---|---|
| 0. Data audit | 1–2 | Frozen cell line table; data audit notebook; medium normalisation table | Go/no-go (Section 4) |
| 1. Network and targets | 2–3 | Human-GEM parsing; metabolite mapping table; gene sets; final target list | Commit `expectations.md` with target list **before** any modelling |
| 2. Pipeline and controls | 4–5 | CV pipeline; permuted-target control; simulated sanity check (a synthetic target built from known genes is recovered) | P control ≈ 0; synthetic target recovered |
| 3. Main analysis | 5–6 | M0/M1/M2, nulls, sensitivity analyses; figures | All leakage checklist items "no" |
| 4. Interpretation and write-up | 7–8 | Case studies; README report; limitations | Results reported for every target |

---

## 14. Deliverables

1. **Figures (core set):**
   - Per-target ΔR²(M1) with the N-met null distribution (violin or interval) and N-all null.
   - ΔR²(M1) vs ΔR²(M2) scatter with capture ratio.
   - Baseline R²(M0) per target, coloured by target class.
   - Negative control (P) distribution.
   - Case-study panels.
2. **Tables:** metabolite mapping table (with exclusion reasons); target gene sets; per-target results.
3. **Write-up:** a concise report (README or short paper-style document) stating the question, design, results for all targets, limitations, and what a negative result means.
4. **Reproducibility:** pinned data versions and Human-GEM commit; environment file; single entry point to rerun analyses from raw files.

Suggested repository layout (to be created during implementation, not now):
```
data/raw/            # downloaded files, not committed
data/processed/      # frozen tables (cell lines, targets, gene sets)
src/                 # parsing, matching, modelling code
notebooks/           # exploration and figures only
results/             # metrics, null distributions, figures
docs/                # these documents and the decision log
```

---

## 15. Decision log

Decisions taken in this plan (to be preserved and extended; each later change gets a dated entry with a reason). **From D10 onward the log is kept in `docs/decision_log.md`.**

| # | Date | Decision | Reason |
|---|---|---|---|
| D1 | 2026-11-16 | Polar metabolites only | Mapping reliability; lipid GPR instability |
| D2 | 2026-11-16 | Solid-tumour lines as primary population | Removes adherent/suspension axis; comparable with Li 2022 |
| D3 | 2026-11-16 | First-neighbour gene sets, transporters included, union of GPR genes | Simplest biologically grounded definition |
| D4 | 2026-11-16 | Target gene-set size 3–60 | Excludes degenerate and hub targets without using outcomes |
| D5 | 2026-11-16 | Ridge with two-stage covariate adjustment | Stable, conservative, interpretable |
| D6 | 2026-11-16 | Random metabolic genes (N-met) as primary null | Separates network locality from general metabolic-gene information |
| D7 | 2026-11-16 | Pin Human-GEM v2.0.0 (commit hash recorded at download) | Most curated release; version drift affects gene sets |
| D8 | 2026-11-16 | CCLE 2019 RNA-seq (RSEM TPM) as expression source | Same release/namespace as metabolomics (alternative: current DepMap release) |
| D9 | 2026-11-16 | No FBA, no multi-omics, no external datasets, no flexible ML in core scope | Keep the project finishable and the question clean |
