# Phase 1 Report: Network, Mapping and Target Selection (FINAL)

**Project:** MetaboNet · **Date:** 2026-11-16 · **Scope:** Phase 1 only (`research_plan.md` §3, §5, §6).
**No model was fitted.** Nothing was computed relating metabolite values to expression. Target selection used only:

- ID mapping;
- the Human-GEM network;
- expression *levels* (the expressed-gene filter);
- each metabolite's own pre-imputation missingness and spread.

**Reproduce:** `make phase1`, then `.venv/bin/python -m src.phase1.report_tables` (regenerates every table below and `expectations.md` §9 from `data/processed/phase1/`).
**Status:** the target list is **frozen** in `expectations.md` §9 after the project owner's review (decisions D31–D38). Phase 2 has not started.

---

## 1. Go/no-go

| Criterion (research_plan §4) | Result | Status |
|---|---|---|
| ≥ 15 polar targets pass rules 1–4 | **74** | **PASS** |

The rules produce more targets than the plan's guess of about 20–35. The number is whatever the pre-specified rules give (research_plan §5). The modelling cost scales linearly; a compute contingency is pre-specified (S5, D38).

## 2. Network: Human-GEM v2.0.0

- **Model:** tag `v2.0.0`, commit `635f533152dc5f7290ce04d12700eaa882273c3e`. The git blob SHA-1 of all six files was checked against the commit (`docs/data_sources.md`).
- **Contents:** 12,931 reactions, 8,461 metabolites (compartment copies), 2,848 genes and 9 compartments, parsed with COBRApy 0.32.1. No solver was used.
- **Excluded reactions (D23, D33):** 1,660 boundary reactions (exchange/sink/demand) and the *Exchange/demand reactions*, *Artificial reactions* (biomass), *Pool reactions* and *Protein degradation* subsystems. Aminoacyl-tRNA synthetase reactions sit in their own subsystem (*Aminoacyl-tRNA biosynthesis*, 21 reactions) and are kept.
- **Annotation tables:** `metabolites.tsv` (KEGG/HMDB/ChEBI per compartment copy), `genes.tsv` (symbols) and `reactions.tsv` were read with pandas.
- **Newer releases:** v2.0.1, v2.1.0 and v2.1.1 exist; v2.0.0 stays pinned (D30).

## 3. Gene ID mapping (GENCODE v19 → Human-GEM)

- **Rule (D26):** Ensembl ID with the version stripped, falling back to the HGNC symbol if it is unique in v19. The 47 `ENSGR` chrY-PAR rows were dropped first (D22).
- **Result:** **2,845 of 2,848 Human-GEM genes mapped (99.9%):** 2,827 by Ensembl ID and 18 by symbol. No CCLE gene is claimed by two Human-GEM genes.
- **Expressed:** 2,123 mapped Human-GEM genes are expressed in the primary set (mean log2(TPM+1) ≥ 1). The matrix as a whole has 15,106 expressed genes.
- **No target loses more than 20% of its genes to unmapped IDs** (research_plan §10.2 flag).
- **Non-coding v19 biotypes:** 15 mapped genes carry one in v19, e.g. *SRD5A2*, *UHRF1* and *ABO* as `processed_transcript`, *ADH1C* and *GBA3* as `polymorphic_pseudogene`. They are kept, since Human-GEM lists them as metabolic genes. The full list is in `gene_mapping.csv`.

**Every unmapped Human-GEM gene:**

{{UNMAPPED_GENES}}

## 4. Metabolite mapping (136 polar features)

**Procedure (D28, refined by D32):**

- **Curation.** One KEGG anchor ID per species (ChEBI or HMDB where KEGG has no entry) in `data/curation/polar_feature_ids.tsv`, with notes on LC-MS isomer issues.
- **Verification.**
  - Every KEGG ID is checked against KEGG's own names. This caught one curation error (cotinine).
  - ChEBI IDs come from KEGG DBLINKS and are expanded across conjugate acid/base forms.
  - HMDB IDs come from ChEBI cross-references (hmdb.ca blocks scripts).
- **Matching.** A feature matches a compartment-free Human-GEM metabolite through any shared KEGG/ChEBI/HMDB ID, **or (D32)** through a name match when the molecular formula agrees and no competing candidate exists.
- **Extra check.** A duplicate scan across Human-GEM found the C18:1-carnitine isomers and the α-NAD anomer.
- **Review.** 23 features triggered review. The project owner's decisions (D31) are recorded in `data/curation/review_decisions.tsv` and applied by the pipeline, so both the rule-based status and the final status are kept.

**Final status:**

{{STATUS_COUNTS}}

- **Ambiguous (19):**
  - the 13 multi-species features;
  - *hexoses (HILIC neg)* and *hexoses (HILIC pos)* (generic hexose: glucose/fructose/galactose/mannose, with α/β-glucose anomers also in Human-GEM);
  - *2-hydroxyglutarate* (D- and L-2HG are both in Human-GEM and are biologically distinct; LC-MS does not resolve them);
  - *oleylcarnitine* (five C18:1-carnitine species in Human-GEM);
  - **lactose** and **pyroglutamic acid**, excluded on review (D31).
- **Unmapped (5), not in Human-GEM:** NMMA, cotinine (a xenobiotic), methionine sulfoxide, heptanoylcarnitine, anserine.
- **Provisional (5), left undecided because they fail rule 3 or 4 regardless:** NAD+, α-hydroxybutyrate, allantoin, lauroylcarnitine, myristoylcarnitine. (Acetylglycine and arachidonylcarnitine, also moot, resolved to confident by the D32 name + formula rule; they still fail rule 3 and rule 4 respectively.)
- **Full mapping table:** Appendix A and `data/processed/phase1/metabolite_mapping.csv`.

### 4.1 Review decisions (D31)

{{REVIEW}}

## 5. Gene sets

Built for every feature with a single Human-GEM metabolite: 112 features (D23, D24, D33).

**Rules:**

1. Collapse compartment copies.
2. Take every reaction involving the metabolite.
3. Drop the excluded reaction classes.
4. Take the union of all GPR genes.

The **enzymatic-only** version also drops reactions where the target itself appears on both sides, i.e. transport. Unit tests are in `tests/test_genesets.py`.

- **Set sizes among the 74 targets:** median 10.5 expressed genes, range 3–59.
- **Rule 3 failures:**
  - 21 features have **< 3** expressed genes, among them 1-methylnicotinamide (2), N-carbamoyl-β-alanine (0), betaine (1), creatinine (1) and thiamine (2);
  - 3 have **> 60**: NADP+ (197), NAD+ (227) and AMP (280). Glutamate (72 before D33) fell to 58 once the *Protein degradation* reactions were excluded and is now a target.
- **Effect of D33 (Protein degradation excluded):** the 18 amino-acid targets lost the 14 shared protease genes (e.g. valine 26 → 12, tryptophan 26 → 12, glycine 51 → 38), and the 15-member amino-acid cluster broke into three small ones (Section 5.1).
- **Other subsystems that put the same genes into > 10 targets:** only *Transport reactions*. Its nine shared genes are genuine first-neighbour transporters (SLC7A5/SLC3A2, SLC7A6/SLC7A7, SLC36A1, SLC38A2, SLC43A1, SLC29A2, SLC7A8) and are kept (D3).

{{SHARED_SUBSYSTEMS}}

- **Reactions without a GPR:** mostly gene-less amino-acid transport reactions. They add no genes, so they only lower the gene count. Per-target counts are in `gene_set_summary.csv`.

### 5.1 Overlap and clusters (D37)

Pairwise Jaccard of the expressed primary sets across the 74 targets: 2,701 pairs, 497 share at least one gene, **33 pairs have Jaccard > 0.5**. Single linkage at Jaccard > 0.5 gives **56 clusters**, five of them with more than one target:

{{CLUSTERS}}

Acetyl-/propionylcarnitine and palmitoyl-/stearoylcarnitine have identical sets (Jaccard = 1.0). The clustered secondary reading of H2 uses these clusters (`target_clusters.csv`); the full-target reading stays primary.

### 5.2 Sensitivity analysis 2 (D36)

Targets with fewer than 3 expressed enzymatic genes are dropped from the enzymatic-only analysis only. **7 of 74** are dropped; 67 are analysed:

{{SENS2_DROPPED}}

## 6. Positive controls (D34, option b)

In Human-GEM, three of the five literature genes are not first neighbours of their metabolite. H5 was amended **before any modelling** to name the in-set genes; M1 and rule 3 are unchanged.

{{POSCONTROLS}}

- **IDO1/TDO2** make *N*-formylkynurenine; AFMID then makes kynurenine (distance 2). **GATM** makes guanidinoacetate; GAMT makes creatine (distance 2). **DPYD** acts one step before DPYS.
- **Controls:** kynurenine (AFMID, KYNU), creatine (GAMT, SLC6A8), taurine (SLC6A6). All three pass selection.
- **Excluded by rule 3:** 1-methylnicotinamide (2 expressed genes: NNMT, AOX1) and N-carbamoyl-β-alanine (0; DPYS and UPB1 in set but not expressed).
- **Pipeline diagnostic, not a target:** M1 for 1-methylnicotinamide with its 2-gene set, with ΔR² and the N-met null, labelled *outside selection rules*; it enters none of H1–H4 (`pipeline_diagnostic_gene_sets.csv`).
- Kynurenine is the concrete motivation for stretch goal 2 (second-neighbour sets).

## 7. Target class (D29, confirmed D35)

Rule: *medium-supplied* = in the RPMI-1640 or DMEM base formulation and not synthesised by human cells; *mixed* = in the medium and synthesised; *synthesised* = in neither. RPMI-1640 and DMEM cover 489 of 735 primary lines. **Exceptions to state in the write-up:** F12-based media (78 lines: DMEM:F12, F12, F10, F12K) supply alanine, aspartate, asparagine, glutamate, proline, putrescine, hypoxanthine and thymidine, which the rule classes as *synthesised*; Leibovitz L-15 (12 lines) replaces glucose with galactose (relevant to sorbitol/galactitol) and supplies alanine and asparagine.

**Among the 74 targets:** 11 medium-supplied, 12 mixed, 51 synthesised.

## 8. Target selection

Rules applied in order (research_plan §5; D18, D25, D27):

1. HILIC (polar) assay;
2. 1:1 confident mapping (after D31/D32);
3. 3–60 expressed genes in the primary gene set;
4. (a) ≤ 20% missing before imputation and (b) SD ≥ 0.05 log10 units.

{{ATTRITION}}

### Frozen target list (identical to `expectations.md` §9)

k = expressed genes in the primary or enzymatic-only gene set. S2 = kept in sensitivity analysis 2. Cluster = multi-target Jaccard cluster (blank = singleton).

{{TARGETS}}

## 9. Outputs (`data/processed/phase1/`, hashes in `FREEZE.tsv`)

| File | Content |
|---|---|
| `metabolite_mapping.csv` | All 136 polar features: species, IDs, Human-GEM candidates and match route, rule-based and final status, reasons, review decision |
| `metabolite_mapping_review_decisions.csv` | The 23 reviewed features: issue, my recommendation, the owner's decision, final status |
| `gene_mapping.csv`, `gene_mapping_unmapped.csv` | Every Human-GEM gene → CCLE gene ID, method, biotype; the 3 unmapped genes |
| `gene_set_summary.csv` | Per feature: reactions (all / excluded / retained / transport), reactions without a GPR, k at each stage, for both versions |
| `gene_set_members.csv` | Every gene in every set: version, symbol, mapped, expressed |
| `target_selection.csv`, `target_attrition.csv` | Rule-by-rule flags, result per feature, sensitivity-2 and diagnostic flags, attrition table |
| `positive_controls.csv`, `pipeline_diagnostic_gene_sets.csv` | Control genes vs first-neighbour sets; the MNA diagnostic set |
| `target_jaccard_matrix.csv`, `target_jaccard_pairs.csv`, `target_clusters.csv` | Overlap between targets; single-linkage clusters |
| `shared_subsystem_genes.csv` | Subsystem/gene pairs contributing a gene to > 10 targets |
| `phase1_summary.json` | Headline numbers |

Curated inputs (committed, hand-written): `data/curation/polar_feature_ids.tsv`, `target_annotations.tsv`, `review_recommendations.tsv`, `review_decisions.tsv`.

## 10. Open items carried into Phase 2

1. **Compute contingency (S5, D38):** time one complete target first; if the projected total exceeds ~8 h, reduce B to 500 for N-met and N-all and log it before running further targets.
2. **HMDB verification:** two HMDB IDs (hexanoyl- and malonylcarnitine) are curator-supplied and unverified (labelled in the table); the mappings rest on name + formula (D32) and, for malonylcarnitine, on Human-GEM's own HMDB annotation.
3. **Sorbitol caveat (D31):** galactitol is possible in the 12 L-15 lines.
4. **Stretch goal 2** (second-neighbour sets) is motivated by kynurenine; it stays a stretch goal.

---

## Appendix A: Mapping table (all 136 polar features)

"Status" shows the final status and, where a review decision applied, the decision. "Selection" is the target-rule outcome for features with a single proposed Human-GEM metabolite.

{{MAPPING}}
