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

| Human-GEM gene | Symbol | Reason |
|---|---|---|
| ENSG00000273841 | TAF9 | symbol TAF9 -> ENSG00000085231.9, which is already the direct Ensembl match of another Human-GEM gene |
| ENSG00000288702 | UGT1A3 | Ensembl ID absent from GENCODE v19 and symbol 'UGT1A3' not a v19 gene_name |
| ENSG00000288705 | UGT1A5 | Ensembl ID absent from GENCODE v19 and symbol 'UGT1A5' not a v19 gene_name |

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

| Status | Features |
|---|---|
| 1:1 confident | 107 |
| 1:1 confident (provisional) | 5 |
| ambiguous | 19 |
| unmapped | 5 |

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

| Feature | Rule-based status | Issue | Decision (D31) | Final status | Rationale | Selection |
|---|---|---|---|---|---|---|
| 3-phosphoglycerate | 1:1 confident (provisional) | LC-MS isomer concern: 2-phosphoglycerate (C00631) is isobaric; not listed as a separate feature | accept | 1:1 confident | 3-PG is the dominant phosphoglycerate pool and Li et al. assigned the feature with a 3-PG standard. 2-PG is isobaric and may co-elute, but PGAM1 (interconverts 3-PG and 2-PG) is already in the gene set, so the set covers the sum of both species. | passes |
| alanine | 1:1 confident (provisional) | enantiomer present in Human-GEM (MAM01638); LC-MS does not resolve enantiomers | accept | 1:1 confident | Enantiomer rule (D31): LC-MS does not resolve L/D; D-amino acids are trace in cultured human cancer cells. Mapped to the L-form (MAM01307); D-alanine (MAM01638) reactions are not added to the gene set. | passes |
| alpha-glycerophosphate | 1:1 confident (provisional) | LC-MS isomer concern: glycerol 2-phosphate (C02979) is isobaric; name implies the alpha (1/3) form | accept | 1:1 confident | 'alpha' names the 1/3 isomer; glycerol 2-phosphate is a minor cellular species. Mapped to sn-glycerol 3-phosphate (MAM02914). | passes |
| arginine | 1:1 confident (provisional) | enantiomer present in Human-GEM (MAM01640); LC-MS does not resolve enantiomers | accept | 1:1 confident | Enantiomer rule (D31). Mapped to L-arginine (MAM01365); D-arginine (MAM01640) reactions not added. | passes |
| aspartate | 1:1 confident (provisional) | enantiomer present in Human-GEM (MAM01641); LC-MS does not resolve enantiomers | accept | 1:1 confident | Enantiomer rule (D31). Mapped to L-aspartate (MAM01370); D-aspartate (MAM01641) reactions not added. | passes |
| hexanoylcarnitine | 1:1 confident | name match MAM03498 (Human-GEM entry has no KEGG/ChEBI/HMDB ID), formula consistent, no competing candidate (D32) | accept | 1:1 confident | Name match to Human-GEM 'Hexanoyl Carnitine' (MAM03498, which carries no external IDs); formula consistent; no competing candidate (D32 rule). | passes |
| lactate | 1:1 confident (provisional) | enantiomer present in Human-GEM (MAM01716); LC-MS does not resolve enantiomers | accept | 1:1 confident | Enantiomer rule (D31). D-lactate (methylglyoxal pathway) is typically <1% of cellular lactate. Mapped to L-lactate (MAM02403); D-lactate (MAM01716) reactions not added. | passes |
| malonylcarnitine | 1:1 confident | single species; unique Human-GEM match MAM02443 via MAM02443(hmdb) | accept | 1:1 confident | Name match to Human-GEM 'malonyl-carnitine' (MAM02443) and shared HMDB0002095; formula consistent; no competing candidate (D32 rule). | passes |
| ornithine | 1:1 confident (provisional) | enantiomer present in Human-GEM (MAM01740); LC-MS does not resolve enantiomers | accept | 1:1 confident | Enantiomer rule (D31). Mapped to L-ornithine (MAM02658); D-ornithine (MAM01740) reactions not added. | passes |
| proline | 1:1 confident (provisional) | enantiomer present in Human-GEM (MAM01742); LC-MS does not resolve enantiomers | accept | 1:1 confident | Enantiomer rule (D31). Mapped to L-proline (MAM02770); D-proline (MAM01742) reactions not added. | passes |
| serine | 1:1 confident (provisional) | enantiomer present in Human-GEM (MAM01744); LC-MS does not resolve enantiomers | accept | 1:1 confident | Enantiomer rule (D31). Mapped to L-serine (MAM02896); D-serine (MAM01744, made by SRR) reactions not added. | passes |
| sorbitol | 1:1 confident (provisional) | LC-MS isomer concern: isobaric with other hexitols (mannitol, galactitol) | accept | 1:1 confident | Sorbitol (AKR1B1/SORD polyol pathway) is the hexitol human cells make; mannitol is not synthesised by human cells. Caveat: the 12 primary lines in Leibovitz L-15 are grown on galactose, where galactitol (AKR1B1 product, isobaric) is possible. | passes |
| stearoylcarnitine | 1:1 confident | name match MAM02940 (Human-GEM entry has no KEGG/ChEBI/HMDB ID), formula consistent, no competing candidate (D32) | accept | 1:1 confident | Name match to Human-GEM 'stearoylcarnitine' (MAM02940, no external IDs); formula consistent; no competing candidate (D32 rule). | passes |
| lactose | 1:1 confident (provisional) | LC-MS isomer concern: isobaric with other disaccharides (e.g. sucrose listed separately) | exclude | ambiguous | Lactose is synthesised only by lactating mammary epithelium (LALBA/B4GALT1); an intracellular signal in cancer lines more plausibly reflects FBS carry-over or an unresolved isobaric disaccharide. A first-neighbour gene set for lactose is not biologically meaningful here. | mapping: ambiguous |
| pyroglutamic acid | 1:1 confident (provisional) | LC-MS isomer concern: can form in-source/ex vivo by cyclisation of glutamine/glutamate | exclude | ambiguous | 5-oxoproline forms non-enzymatically and in the ion source by cyclisation of glutamine and glutamate, so the measured signal is partly an artefact of much larger pools. | mapping: ambiguous |
| oleylcarnitine | ambiguous | Human-GEM has five C18:1 acylcarnitines: L-oleoylcarnitine (MAM02410, no IDs), octadecenoylcarnitine(5) (MAM02639, HMDB0005065), (7) (MAM02640), (11)  | keep_ambiguous | ambiguous | Five C18:1-carnitine species in Human-GEM; the LC-MS transition does not distinguish double-bond position. | mapping: ambiguous |
| NAD | 1:1 confident (provisional) | Human-GEM also contains alpha-NAD(+) (MAM20063), an anomer of NAD+ (MAM02552) | moot | 1:1 confident (provisional) | Not decided: fails rule 3 (227 expressed genes) regardless. alpha-NAD+ (MAM20063) is a minor anomer. | gene set 227 expressed genes (> 60) |
| acetylglycine | 1:1 confident | name match MAM03397 (Human-GEM entry has no KEGG/ChEBI/HMDB ID), formula consistent, no competing candidate (D32) | moot | 1:1 confident | Not decided: fails rule 3 (0 expressed genes) regardless. | gene set 0 expressed genes (< 3) |
| allantoin | 1:1 confident (provisional) | name-only match MAM01313; formula consistent; curated (S)-allantoin (C02350); Human-GEM allantoin is the unspecified/racemic form (C01551, CHEBI:15676 | moot | 1:1 confident (provisional) | Not decided: fails rule 3 (0 expressed genes) regardless; humans lack uricase. | gene set 0 expressed genes (< 3) |
| alpha-hydroxybutyrate | 1:1 confident (provisional) | LC-MS isomer concern: LC-MS does not resolve R/S; isobaric with 3-hydroxybutyrate (separation not documented) | moot | 1:1 confident (provisional) | Not decided: fails rule 4 (37% missing) regardless. | 37% missing |
| arachidonyl_carnitine | 1:1 confident | name match MAM01363 (Human-GEM entry has no KEGG/ChEBI/HMDB ID), formula consistent, no competing candidate (D32) | moot | 1:1 confident | Not decided: fails rule 4 (32% missing) regardless. | 32% missing |
| lauroylcarnitine | 1:1 confident (provisional) | match ['MAM01729'] rests only on a curator-supplied HMDB ID that could not be verified | moot | 1:1 confident (provisional) | Not decided: fails rule 4 (22% missing) regardless. | 22% missing |
| myristoylcarnitine | 1:1 confident (provisional) | match ['MAM02973'] rests only on a curator-supplied HMDB ID that could not be verified | moot | 1:1 confident (provisional) | Not decided: fails rule 4 (21% missing) regardless. | 21% missing |

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

| Subsystem | Genes in > 10 targets | Max targets | Genes (n targets) |
|---|---|---|---|
| Transport reactions | 9 | 20 | SLC36A1 (20), SLC7A5 (19), SLC7A6 (19), SLC3A2 (18), SLC7A7 (17), SLC7A8 (16), SLC43A1 (14), SLC29A2 (11), SLC38A2 (11) |

- **Reactions without a GPR:** mostly gene-less amino-acid transport reactions. They add no genes, so they only lower the gene count. Per-target counts are in `gene_set_summary.csv`.

### 5.1 Overlap and clusters (D37)

Pairwise Jaccard of the expressed primary sets across the 74 targets: 2,701 pairs, 497 share at least one gene, **33 pairs have Jaccard > 0.5**. Single linkage at Jaccard > 0.5 gives **56 clusters**, five of them with more than one target:

| Cluster | Targets (n) | Members |
|---|---|---|
| C01 | 9 | 2-deoxyadenosine, 2-deoxycytidine, adenosine, cytidine, guanosine, inosine, thymidine, uridine, xanthosine |
| C02 | 6 | acetylcarnitine, carnitine, malonylcarnitine, palmitoylcarnitine, propionylcarnitine, stearoylcarnitine |
| C03 | 3 | isoleucine, leucine, valine |
| C04 | 3 | phenylalanine, tryptophan, tyrosine |
| C05 | 2 | asparagine, threonine |

Acetyl-/propionylcarnitine and palmitoyl-/stearoylcarnitine have identical sets (Jaccard = 1.0). The clustered secondary reading of H2 uses these clusters (`target_clusters.csv`); the full-target reading stays primary.

### 5.2 Sensitivity analysis 2 (D36)

Targets with fewer than 3 expressed enzymatic genes are dropped from the enzymatic-only analysis only. **7 of 74** are dropped; 67 are analysed:

| Target | k (primary) | k (enzymatic) |
|---|---|---|
| glucuronate | 3 | 2 |
| oxalate | 6 | 1 |
| thymine | 3 | 2 |
| histidine | 14 | 2 |
| tryptophan | 12 | 2 |
| taurine | 3 | 1 |
| carnosine | 4 | 2 |

## 6. Positive controls (D34, option b)

In Human-GEM, three of the five literature genes are not first neighbours of their metabolite. H5 was amended **before any modelling** to name the in-set genes; M1 and rule 3 are unchanged.

| Control | Literature gene(s) | H5 in-set gene(s) | k expressed | Role |
|---|---|---|---|---|
| 1-methylnicotinamide | NNMT:expressed | NNMT:expressed | 2 | pipeline diagnostic (not a target) |
| kynurenine | IDO1:NOT in set, TDO2:NOT in set | AFMID:expressed, KYNU:expressed | 7 | target + H5 control |
| creatine | SLC6A8:expressed, GATM:NOT in set | GAMT:expressed, SLC6A8:expressed | 6 | target + H5 control |
| N-carbamoyl-beta-alanine | DPYD:NOT in set | DPYS:in set, not expressed, UPB1:in set, not expressed | 0 | excluded by rule 3 |
| taurine | SLC6A6:expressed | SLC6A6:expressed | 3 | target + H5 control |

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

| Step | n |
|---|---|
| CCLE metabolomics features (all assays) | 225 |
| Rule 1: measured by HILIC (polar) | 136 |
| Rule 2: 1:1 confident mapping | 107 |
| Rule 3: 3-60 expressed genes (primary gene set) | 85 |
| Rule 4a: <= 20% missing before imputation | 74 |
| Rule 4b: SD >= 0.05 log10 units | 74 |

### Frozen target list (identical to `expectations.md` §9)

k = expressed genes in the primary or enzymatic-only gene set. S2 = kept in sensitivity analysis 2. Cluster = multi-target Jaccard cluster (blank = singleton).

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

| Feature | Status | KEGG | Human-GEM | Reason | Selection |
|---|---|---|---|---|---|
| 2-aminoadipate | 1:1 confident | C00956 | MAM02322 | single species; unique Human-GEM match MAM02322 via MAM02322(kegg) | gene set 2 expressed genes (< 3) |
| 3-phosphoglycerate | 1:1 confident (reviewed: accept) | C00197 | MAM00913 | 3-PG is the dominant phosphoglycerate pool and Li et al. assigned the feature with a 3-PG standard. 2-PG is isobaric and may co-elute, but PGAM1 (interconverts 3-PG and 2-PG) is al | passes |
| alpha-glycerophosphate | 1:1 confident (reviewed: accept) | C00093 | MAM02914 | 'alpha' names the 1/3 isomer; glycerol 2-phosphate is a minor cellular species. Mapped to sn-glycerol 3-phosphate (MAM02914). | passes |
| 4-pyridoxate | 1:1 confident | C00847 | MAM01033 | single species; unique Human-GEM match MAM01033 via MAM01033(chebi/hmdb/kegg) | gene set 1 expressed genes (< 3); 29% missing |
| aconitate | 1:1 confident | C00417 | MAM01580 | single species; unique Human-GEM match MAM01580 via MAM01580(chebi/kegg) | gene set 2 expressed genes (< 3) |
| adenine | 1:1 confident | C00147 | MAM01279 | single species; unique Human-GEM match MAM01279 via MAM01279(chebi/hmdb/kegg) | 23% missing |
| adipate | 1:1 confident | C06104 | MAM03408 | single species; unique Human-GEM match MAM03408 via MAM03408(chebi/hmdb/kegg) | gene set 1 expressed genes (< 3) |
| alpha-ketoglutarate | 1:1 confident | C00026 | MAM01306 | single species; unique Human-GEM match MAM01306 via MAM01306(chebi/hmdb/kegg) | passes |
| AMP | 1:1 confident | C00020 | MAM01334 | single species; unique Human-GEM match MAM01334 via MAM01334(chebi/hmdb/kegg) | gene set 280 expressed genes (> 60) |
| citrate | 1:1 confident | C00158 | MAM01587 | single species; unique Human-GEM match MAM01587 via MAM01587(chebi/hmdb/kegg) | passes |
| isocitrate | 1:1 confident | C00311 | MAM02183 | single species; unique Human-GEM match MAM02183 via MAM02183(chebi/hmdb/kegg) | passes |
| CMP | 1:1 confident | C00055 | MAM01590 | single species; unique Human-GEM match MAM01590 via MAM01590(chebi/hmdb/kegg) | passes |
| cystathionine | 1:1 confident | C02291 | MAM02349 | single species; unique Human-GEM match MAM02349 via MAM02349(chebi/hmdb/kegg) | passes |
| cytidine | 1:1 confident | C00475 | MAM01630 | single species; unique Human-GEM match MAM01630 via MAM01630(chebi/hmdb/kegg) | passes |
| dCMP | 1:1 confident | C00239 | MAM01644 | single species; unique Human-GEM match MAM01644 via MAM01644(chebi/kegg) | 30% missing |
| DHAP/glyceraldehyde 3P | ambiguous | C00111; C00118 | MAM01690;MAM01939 | summed/multi-species LC-MS feature (2 species; glycerone phosphate: 1 HG match(es), D-glyceraldehyde 3-phosphate: 1 HG match(es)) | mapping: ambiguous |
| erythrose-4-phosphate | 1:1 confident | C00279 | MAM01785 | single species; unique Human-GEM match MAM01785 via MAM01785(chebi/kegg) | passes |
| F1P/F6P/G1P/G6P | ambiguous | C01094; C00085; C00103; C00092 | MAM01842;MAM01845;MAM01967;MAM01968 | summed/multi-species LC-MS feature (4 species; D-fructose 1-phosphate: 1 HG match(es), D-fructose 6-phosphate: 1 HG match(es), D-glucose 1-phosphate: 1 HG match(es), D-glucose 6-ph | mapping: ambiguous |
| hexoses (HILIC neg) | ambiguous | C00031; C00095; C00124; C00159 | MAM01840;MAM01965;MAM02453 | summed/multi-species LC-MS feature (4 species; D-glucose: 1 HG match(es), D-fructose: 1 HG match(es), D-galactose: 0 HG match(es), D-mannose: 1 HG match(es)) | mapping: ambiguous |
| fumarate/maleate/alpha-ketoisovalerate | ambiguous | C00122; C01384; C00141 | MAM00824;MAM01862 | summed/multi-species LC-MS feature (3 species; fumarate: 1 HG match(es), maleate: 0 HG match(es), 3-methyl-2-oxobutanoate: 1 HG match(es)) | mapping: ambiguous |
| glucuronate | 1:1 confident | C00191 | MAM01973 | single species; unique Human-GEM match MAM01973 via MAM01973(hmdb/kegg) | passes |
| glutathione oxidized | 1:1 confident | C00127 | MAM02027 | single species; unique Human-GEM match MAM02027 via MAM02027(chebi/hmdb/kegg) | passes |
| glutathione reduced | 1:1 confident | C00051 | MAM02026 | single species; unique Human-GEM match MAM02026 via MAM02026(chebi/hmdb/kegg) | passes |
| GMP | 1:1 confident | C00144 | MAM02016 | single species; unique Human-GEM match MAM02016 via MAM02016(chebi/kegg) | passes |
| guanosine | 1:1 confident | C00387 | MAM02038 | single species; unique Human-GEM match MAM02038 via MAM02038(chebi/hmdb/kegg) | passes |
| hippurate | 1:1 confident | C01586 | MAM02123 | single species; unique Human-GEM match MAM02123 via MAM02123(chebi/hmdb/kegg) | gene set 2 expressed genes (< 3) |
| hypoxanthine | 1:1 confident | C00262 | MAM02159 | single species; unique Human-GEM match MAM02159 via MAM02159(chebi/hmdb/kegg) | passes |
| inosine | 1:1 confident | C00294 | MAM02170 | single species; unique Human-GEM match MAM02170 via MAM02170(chebi/hmdb/kegg) | passes |
| kynurenine | 1:1 confident | C00328 | MAM02319 | single species; unique Human-GEM match MAM02319 via MAM02319(chebi/hmdb/kegg) | passes |
| lactate | 1:1 confident (reviewed: accept) | C00186 | MAM02403 | Enantiomer rule (D31). D-lactate (methylglyoxal pathway) is typically <1% of cellular lactate. Mapped to L-lactate (MAM02403); D-lactate (MAM01716) reactions not added. | passes |
| lactose | ambiguous (reviewed: exclude) | C00243 | MAM02332 | Lactose is synthesised only by lactating mammary epithelium (LALBA/B4GALT1); an intracellular signal in cancer lines more plausibly reflects FBS carry-over or an unresolved isobari | mapping: ambiguous |
| malate | 1:1 confident | C00149 | MAM02439 | single species; unique Human-GEM match MAM02439 via MAM02439(chebi/hmdb/kegg) | passes |
| NAD | 1:1 confident (provisional) (reviewed: moot) | C00003 | MAM02552 | Not decided: fails rule 3 (227 expressed genes) regardless. alpha-NAD+ (MAM20063) is a minor anomer. | gene set 227 expressed genes (> 60) |
| NADP | 1:1 confident | C00006 | MAM02554 | single species; unique Human-GEM match MAM02554 via MAM02554(chebi/kegg) | gene set 197 expressed genes (> 60) |
| oxalate | 1:1 confident | C00209 | MAM02661 | single species; unique Human-GEM match MAM02661 via MAM02661(chebi/kegg) | passes |
| pantothenate | 1:1 confident | C00864 | MAM02680 | single species; unique Human-GEM match MAM02680 via MAM02680(chebi/hmdb/kegg) | passes |
| PEP | 1:1 confident | C00074 | MAM02696 | single species; unique Human-GEM match MAM02696 via MAM02696(chebi/hmdb/kegg) | passes |
| ribose-5-P/ribulose5-P | ambiguous | C00117; C00199 | MAM02845;MAM02846 | summed/multi-species LC-MS feature (2 species; D-ribose 5-phosphate: 1 HG match(es), D-ribulose 5-phosphate: 1 HG match(es)) | mapping: ambiguous |
| sorbitol | 1:1 confident (reviewed: accept) | C00794 | MAM01682 | Sorbitol (AKR1B1/SORD polyol pathway) is the hexitol human cells make; mannitol is not synthesised by human cells. Caveat: the 12 primary lines in Leibovitz L-15 are grown on galac | passes |
| succinate/methylmalonate | ambiguous | C00042; C02170 | MAM02479;MAM02943 | summed/multi-species LC-MS feature (2 species; succinate: 1 HG match(es), methylmalonate: 1 HG match(es)) | mapping: ambiguous |
| sucrose | 1:1 confident | C00089 | MAM02945 | single species; unique Human-GEM match MAM02945 via MAM02945(chebi/hmdb/kegg) | gene set 1 expressed genes (< 3) |
| thymine | 1:1 confident | C00178 | MAM02997 | single species; unique Human-GEM match MAM02997 via MAM02997(chebi/kegg) | passes |
| UMP | 1:1 confident | C00105 | MAM03114 | single species; unique Human-GEM match MAM03114 via MAM03114(chebi/kegg) | 22% missing |
| UDP-galactose/UDP-glucose | ambiguous | C00052; C00029 | MAM03107;MAM03108 | summed/multi-species LC-MS feature (2 species; UDP-alpha-D-galactose: 1 HG match(es), UDP-alpha-D-glucose: 1 HG match(es)) | mapping: ambiguous |
| uracil | 1:1 confident | C00106 | MAM03118 | single species; unique Human-GEM match MAM03118 via MAM03118(chebi/hmdb/kegg) | passes |
| urate | 1:1 confident | C00366 | MAM03120 | single species; unique Human-GEM match MAM03120 via MAM03120(chebi/hmdb/kegg) | passes |
| uridine | 1:1 confident | C00299 | MAM03123 | single species; unique Human-GEM match MAM03123 via MAM03123(chebi/hmdb/kegg) | passes |
| xanthine | 1:1 confident | C00385 | MAM03148 | single species; unique Human-GEM match MAM03148 via MAM03148(chebi/hmdb/kegg) | passes |
| taurocholate | 1:1 confident | C05122 | MAM02963 | single species; unique Human-GEM match MAM02963 via MAM02963(chebi/hmdb/kegg) | 36% missing |
| glycodeoxycholate/glycochenodeoxycholate | ambiguous | C05464; C05466 | MAM01987;MAM01989 | summed/multi-species LC-MS feature (2 species; glycodeoxycholate: 1 HG match(es), glycochenodeoxycholate: 1 HG match(es)) | mapping: ambiguous |
| taurodeoxycholate/taurochenodeoxycholate | ambiguous | C05463; C05465 | MAM02962;MAM02964 | summed/multi-species LC-MS feature (2 species; taurodeoxycholate: 1 HG match(es), taurochenodeoxycholate: 1 HG match(es)) | mapping: ambiguous |
| phosphocreatine | 1:1 confident | C02305 | MAM01620 | single species; unique Human-GEM match MAM01620 via MAM01620(chebi/kegg) | passes |
| 3-methyladipate/pimelate | ambiguous | ; C02656 | — | summed/multi-species LC-MS feature (2 species; 3-methyladipate: 0 HG match(es), pimelate: 0 HG match(es)) | mapping: ambiguous |
| 6-phosphogluconate | 1:1 confident | C00345 | MAM01169 | single species; unique Human-GEM match MAM01169 via MAM01169(chebi/hmdb/kegg) | 24% missing |
| alpha-hydroxybutyrate | 1:1 confident (provisional) (reviewed: moot) | C05984 | MAM00648 | Not decided: fails rule 4 (37% missing) regardless. | 37% missing |
| 2-hydroxyglutarate | ambiguous | C03196; C01087 | MAM20012 | summed/multi-species LC-MS feature (2 species; (S)-2-hydroxyglutarate: 0 HG match(es), (R)-2-hydroxyglutarate: 1 HG match(es)) | mapping: ambiguous |
| inositol | 1:1 confident | C00137 | MAM02171 | single species; unique Human-GEM match MAM02171 via MAM02171(chebi/hmdb/kegg) | passes |
| malondialdehyde | 1:1 confident | C19440 | MAM02441 | single species; unique Human-GEM match MAM02441 via MAM02441(chebi/kegg) | gene set 1 expressed genes (< 3) |
| glycine | 1:1 confident | C00037 | MAM01986 | single species; unique Human-GEM match MAM01986 via MAM01986(chebi/hmdb/kegg) | passes |
| alanine | 1:1 confident (reviewed: accept) | C00041 | MAM01307 | Enantiomer rule (D31): LC-MS does not resolve L/D; D-amino acids are trace in cultured human cancer cells. Mapped to the L-form (MAM01307); D-alanine (MAM01638) reactions are not a | passes |
| serine | 1:1 confident (reviewed: accept) | C00065 | MAM02896 | Enantiomer rule (D31). Mapped to L-serine (MAM02896); D-serine (MAM01744, made by SRR) reactions not added. | passes |
| threonine | 1:1 confident | C00188 | MAM02993 | single species; unique Human-GEM match MAM02993 via MAM02993(chebi/hmdb/kegg) | passes |
| methionine | 1:1 confident | C00073 | MAM02471 | single species; unique Human-GEM match MAM02471 via MAM02471(chebi/hmdb/kegg) | passes |
| aspartate | 1:1 confident (reviewed: accept) | C00049 | MAM01370 | Enantiomer rule (D31). Mapped to L-aspartate (MAM01370); D-aspartate (MAM01641) reactions not added. | passes |
| glutamate | 1:1 confident | C00025 | MAM01974 | single species; unique Human-GEM match MAM01974 via MAM01974(chebi/hmdb/kegg) | passes |
| asparagine | 1:1 confident | C00152 | MAM01369 | single species; unique Human-GEM match MAM01369 via MAM01369(chebi/hmdb/kegg) | passes |
| glutamine | 1:1 confident | C00064 | MAM01975 | single species; unique Human-GEM match MAM01975 via MAM01975(chebi/hmdb/kegg) | passes |
| histidine | 1:1 confident | C00135 | MAM02125 | single species; unique Human-GEM match MAM02125 via MAM02125(chebi/hmdb/kegg) | passes |
| arginine | 1:1 confident (reviewed: accept) | C00062 | MAM01365 | Enantiomer rule (D31). Mapped to L-arginine (MAM01365); D-arginine (MAM01640) reactions not added. | passes |
| lysine | 1:1 confident | C00047 | MAM02426 | single species; unique Human-GEM match MAM02426 via MAM02426(chebi/hmdb/kegg) | passes |
| valine | 1:1 confident | C00183 | MAM03135 | single species; unique Human-GEM match MAM03135 via MAM03135(chebi/hmdb/kegg) | passes |
| leucine | 1:1 confident | C00123 | MAM02360 | single species; unique Human-GEM match MAM02360 via MAM02360(chebi/hmdb/kegg) | passes |
| isoleucine | 1:1 confident | C00407 | MAM02184 | single species; unique Human-GEM match MAM02184 via MAM02184(chebi/hmdb/kegg) | passes |
| phenylalanine | 1:1 confident | C00079 | MAM02724 | single species; unique Human-GEM match MAM02724 via MAM02724(chebi/hmdb/kegg) | passes |
| tyrosine | 1:1 confident | C00082 | MAM03101 | single species; unique Human-GEM match MAM03101 via MAM03101(chebi/hmdb/kegg) | passes |
| tryptophan | 1:1 confident | C00078 | MAM03089 | single species; unique Human-GEM match MAM03089 via MAM03089(chebi/hmdb/kegg) | passes |
| proline | 1:1 confident (reviewed: accept) | C00148 | MAM02770 | Enantiomer rule (D31). Mapped to L-proline (MAM02770); D-proline (MAM01742) reactions not added. | passes |
| cis/trans-hydroxyproline | ambiguous | C01157; C01015 | MAM03037 | summed/multi-species LC-MS feature (2 species; trans-4-hydroxy-L-proline: 1 HG match(es), cis-4-hydroxy-L-proline: 0 HG match(es)) | mapping: ambiguous |
| ornithine | 1:1 confident (reviewed: accept) | C00077 | MAM02658 | Enantiomer rule (D31). Mapped to L-ornithine (MAM02658); D-ornithine (MAM01740) reactions not added. | passes |
| citrulline | 1:1 confident | C00327 | MAM01588 | single species; unique Human-GEM match MAM01588 via MAM01588(chebi/hmdb/kegg) | passes |
| taurine | 1:1 confident | C00245 | MAM02961 | single species; unique Human-GEM match MAM02961 via MAM02961(chebi/hmdb/kegg) | passes |
| 5-HIAA | 1:1 confident | C05635 | MAM01103 | single species; unique Human-GEM match MAM01103 via MAM01103(chebi/hmdb/kegg) | 25% missing |
| serotonin | 1:1 confident | C00780 | MAM02897 | single species; unique Human-GEM match MAM02897 via MAM02897(chebi/hmdb/kegg) | 22% missing |
| GABA | 1:1 confident | C00334 | MAM00970 | single species; unique Human-GEM match MAM00970 via MAM00970(chebi/hmdb/kegg) | passes |
| acetylglycine | 1:1 confident (reviewed: moot) | — | MAM03397 | Not decided: fails rule 3 (0 expressed genes) regardless. | gene set 0 expressed genes (< 3) |
| dimethylglycine | 1:1 confident | C01026 | MAM01708 | single species; unique Human-GEM match MAM01708 via MAM01708(chebi/hmdb/kegg) | gene set 1 expressed genes (< 3) |
| homocysteine | 1:1 confident | C00155 | MAM02133 | single species; unique Human-GEM match MAM02133 via MAM02133(hmdb/kegg) | passes |
| SDMA/ADMA | ambiguous | ; C03626 | MAM02577 | summed/multi-species LC-MS feature (2 species; symmetric dimethylarginine: 0 HG match(es), asymmetric dimethylarginine: 1 HG match(es)) | mapping: ambiguous |
| NMMA | unmapped | C03884 | — | no Human-GEM metabolite shares a KEGG/ChEBI/HMDB ID or name | mapping: unmapped |
| allantoin | 1:1 confident (provisional) (reviewed: moot) | C02350 | MAM01313 | Not decided: fails rule 3 (0 expressed genes) regardless; humans lack uricase. | gene set 0 expressed genes (< 3) |
| anthranilic acid | 1:1 confident | C00108 | MAM01342 | single species; unique Human-GEM match MAM01342 via MAM01342(chebi/hmdb/kegg) | gene set 2 expressed genes (< 3) |
| kynurenic acid | 1:1 confident | C01717 | MAM00990 | single species; unique Human-GEM match MAM00990 via MAM00990(chebi/hmdb/kegg) | passes |
| 5-adenosylhomocysteine | 1:1 confident | C00021 | MAM02871 | single species; unique Human-GEM match MAM02871 via MAM02871(chebi/hmdb/kegg) | passes |
| carnosine | 1:1 confident | C00386 | MAM01423 | single species; unique Human-GEM match MAM01423 via MAM01423(chebi/hmdb/kegg) | passes |
| N-carbamoyl-beta-alanine | 1:1 confident | C02642 | MAM00923 | single species; unique Human-GEM match MAM00923 via MAM00923(chebi/hmdb/kegg) | gene set 0 expressed genes (< 3) |
| thiamine | 1:1 confident | C00378 | MAM02982 | single species; unique Human-GEM match MAM02982 via MAM02982(chebi/hmdb/kegg) | gene set 2 expressed genes (< 3) |
| niacinamide | 1:1 confident | C00153 | MAM02583 | single species; unique Human-GEM match MAM02583 via MAM02583(chebi/kegg) | passes |
| betaine | 1:1 confident | C00719 | MAM01393 | single species; unique Human-GEM match MAM01393 via MAM01393(chebi/hmdb/kegg) | gene set 1 expressed genes (< 3) |
| choline | 1:1 confident | C00114 | MAM01513 | single species; unique Human-GEM match MAM01513 via MAM01513(chebi/hmdb/kegg) | passes |
| alpha-glycerophosphocholine | 1:1 confident | C00670 | MAM02912 | single species; unique Human-GEM match MAM02912 via MAM02912(chebi/hmdb/kegg) | passes |
| acetylcholine | 1:1 confident | C01996 | MAM01260 | single species; unique Human-GEM match MAM01260 via MAM01260(chebi/kegg) | gene set 2 expressed genes (< 3); 36% missing |
| creatine | 1:1 confident | C00300 | MAM01619 | single species; unique Human-GEM match MAM01619 via MAM01619(chebi/hmdb/kegg) | passes |
| creatinine | 1:1 confident | C00791 | MAM01621 | single species; unique Human-GEM match MAM01621 via MAM01621(chebi/hmdb/kegg) | gene set 1 expressed genes (< 3) |
| thyroxine | 1:1 confident | C01829 | MAM02998 | single species; unique Human-GEM match MAM02998 via MAM02998(chebi/kegg) | 22% missing |
| trimethylamine-N-oxide | 1:1 confident | C01104 | MAM03054 | single species; unique Human-GEM match MAM03054 via MAM03054(chebi/kegg) | gene set 2 expressed genes (< 3) |
| hexoses (HILIC pos) | ambiguous | C00031; C00095; C00124; C00159 | MAM01840;MAM01965;MAM02453 | summed/multi-species LC-MS feature (4 species; D-glucose: 1 HG match(es), D-fructose: 1 HG match(es), D-galactose: 0 HG match(es), D-mannose: 1 HG match(es)) | mapping: ambiguous |
| adenosine | 1:1 confident | C00212 | MAM01280 | single species; unique Human-GEM match MAM01280 via MAM01280(chebi/hmdb/kegg) | passes |
| thymidine | 1:1 confident | C00214 | MAM02996 | single species; unique Human-GEM match MAM02996 via MAM02996(chebi/hmdb/kegg) | passes |
| xanthosine | 1:1 confident | C01762 | MAM03149 | single species; unique Human-GEM match MAM03149 via MAM03149(chebi/hmdb/kegg) | passes |
| 2-deoxyadenosine | 1:1 confident | C00559 | MAM01666 | single species; unique Human-GEM match MAM01666 via MAM01666(chebi/hmdb/kegg) | passes |
| 2-deoxycytidine | 1:1 confident | C00881 | MAM01668 | single species; unique Human-GEM match MAM01668 via MAM01668(chebi/hmdb/kegg) | passes |
| cAMP | 1:1 confident | C00575 | MAM01419 | single species; unique Human-GEM match MAM01419 via MAM01419(chebi/hmdb/kegg) | 35% missing |
| cotinine | unmapped | — | — | no Human-GEM metabolite shares a KEGG/ChEBI/HMDB ID or name | mapping: unmapped |
| pipecolic acid | 1:1 confident | C00408 | MAM02413 | single species; unique Human-GEM match MAM02413 via MAM02413(chebi/hmdb/kegg) | gene set 0 expressed genes (< 3) |
| pyroglutamic acid | ambiguous (reviewed: exclude) | C01879 | MAM01127 | 5-oxoproline forms non-enzymatically and in the ion source by cyclisation of glutamine and glutamate, so the measured signal is partly an artefact of much larger pools. | mapping: ambiguous |
| 1-methylnicotinamide | 1:1 confident | C02918 | MAM00536 | single species; unique Human-GEM match MAM00536 via MAM00536(chebi/hmdb/kegg) | gene set 2 expressed genes (< 3) |
| butyrobetaine | 1:1 confident | C01181 | MAM01922 | single species; unique Human-GEM match MAM01922 via MAM01922(chebi/kegg) | gene set 1 expressed genes (< 3) |
| putrescine | 1:1 confident | C00134 | MAM02812 | single species; unique Human-GEM match MAM02812 via MAM02812(chebi/hmdb/kegg) | passes |
| methionine sulfoxide | unmapped | C02989 | — | no Human-GEM metabolite shares a KEGG/ChEBI/HMDB ID or name | mapping: unmapped |
| carnitine | 1:1 confident | C00318 | MAM02348 | single species; unique Human-GEM match MAM02348 via MAM02348(hmdb) | passes |
| acetylcarnitine | 1:1 confident | C02571 | MAM02634 | single species; unique Human-GEM match MAM02634 via MAM02634(chebi/hmdb/kegg) | passes |
| propionylcarnitine | 1:1 confident | C03017 | MAM02657 | single species; unique Human-GEM match MAM02657 via MAM02657(chebi/hmdb/kegg) | passes |
| malonylcarnitine | 1:1 confident (reviewed: accept) | — | MAM02443 | Name match to Human-GEM 'malonyl-carnitine' (MAM02443) and shared HMDB0002095; formula consistent; no competing candidate (D32 rule). | passes |
| butyrylcarnitine/isobutyrylcarnitine | ambiguous | ;  | MAM02635 | summed/multi-species LC-MS feature (2 species; O-butanoyl-L-carnitine: 1 HG match(es), O-isobutyryl-L-carnitine: 0 HG match(es)) | mapping: ambiguous |
| valerylcarnitine/isovalerylcarnitine/2-methylbutyroylcarnitine | ambiguous | ; C20826;  | — | summed/multi-species LC-MS feature (3 species; O-valeryl-L-carnitine: 0 HG match(es), O-isovaleryl-L-carnitine: 0 HG match(es), O-2-methylbutyryl-L-carnitine: 0 HG match(es)) | mapping: ambiguous |
| hexanoylcarnitine | 1:1 confident (reviewed: accept) | — | MAM03498 | Name match to Human-GEM 'Hexanoyl Carnitine' (MAM03498, which carries no external IDs); formula consistent; no competing candidate (D32 rule). | passes |
| heptanoylcarnitine | unmapped | — | — | no Human-GEM metabolite shares a KEGG/ChEBI/HMDB ID or name | mapping: unmapped |
| lauroylcarnitine | 1:1 confident (provisional) (reviewed: moot) | — | MAM01729 | Not decided: fails rule 4 (22% missing) regardless. | 22% missing |
| myristoylcarnitine | 1:1 confident (provisional) (reviewed: moot) | — | MAM02973 | Not decided: fails rule 4 (21% missing) regardless. | 21% missing |
| palmitoylcarnitine | 1:1 confident | C02990 | MAM02411 | single species; unique Human-GEM match MAM02411 via MAM02411(chebi/hmdb/kegg) | passes |
| stearoylcarnitine | 1:1 confident (reviewed: accept) | — | MAM02940 | Name match to Human-GEM 'stearoylcarnitine' (MAM02940, no external IDs); formula consistent; no competing candidate (D32 rule). | passes |
| oleylcarnitine | ambiguous (reviewed: keep_ambiguous) | — | MAM02639 | Five C18:1-carnitine species in Human-GEM; the LC-MS transition does not distinguish double-bond position. | mapping: ambiguous |
| arachidonyl_carnitine | 1:1 confident (reviewed: moot) | — | MAM01363 | Not decided: fails rule 4 (32% missing) regardless. | 32% missing |
| sarcosine | 1:1 confident | C00213 | MAM02880 | single species; unique Human-GEM match MAM02880 via MAM02880(chebi/hmdb/kegg) | gene set 2 expressed genes (< 3); 35% missing |
| beta-alanine | 1:1 confident | C00099 | MAM01383 | single species; unique Human-GEM match MAM01383 via MAM01383(chebi/hmdb/kegg) | 22% missing |
| anserine | unmapped | C01262 | — | no Human-GEM metabolite shares a KEGG/ChEBI/HMDB ID or name | mapping: unmapped |
