# Literature Review

**Project:** Predicting metabolite abundance from gene expression using human metabolic network information (CCLE)
**Status:** Draft v0.1 (16 November 2026)
**Scope of this document:** What prior work established, where studies disagree, how good the data and tools are, and what is still worth testing. Methodological decisions are in `research_plan.md`. Testable hypotheses are in `expectations.md`.

---

## 1. How this review was assembled

- Primary sources were read directly where accessible (full text for Li et al. 2019, Li et al. 2022 preprint, Cavicchioli et al. 2022, Benedetti et al. 2023, the Recon8D preprint, Cherkaoui et al. 2022, Luo et al. 2026; abstract or summary level for Peng et al. 2018, Ortmayr et al. 2019, Zelezniak et al. 2014, Venet et al. 2011, Vande Voorde et al. 2019 and Joly et al. 2021).
- **Naming note.** Two different first authors share the surname Li. "Li et al. (2019)" is Haoxin Li et al., *Nature Medicine* (the CCLE metabolomics paper). "Li et al. (2022)" is taken here to mean Huaping Li et al., *FASEB Journal*, "Gene expression is a poor predictor of steady-state metabolite abundance in cancer cells", because it is the 2022 study that addresses this exact question with CCLE data. If a different 2022 paper was intended, Section 3.1 should be revisited.
- **Version note.** The published version of the Recon8D work (renamed Metab8D, *Communications Biology*, June 2026) could not be retrieved. Numbers quoted below come from the 2024 bioRxiv preprint and may differ in the final article.
- Where a paper's internal numbers are inconsistent, this is stated rather than silently corrected.

---

## 2. The dataset: CCLE metabolomics (Li et al., 2019)

**What was measured.** 225 metabolites in 928 cancer cell lines from more than 20 cancer types, by LC-MS. The panel comprises **124 polar metabolites** (two HILIC methods: acidic/positive and basic/negative ion mode) and **101 lipids** (reversed-phase LC), measured on triple-quadrupole instruments. This is a **targeted** panel, not untargeted discovery metabolomics.

**Quality control.** Samples were run in four batches. A pooled reference made from 11 cell lines was injected after every 20 samples for drift and batch correction, followed by median normalisation. Identity was checked against known standards and manual peak review; formal identification-confidence levels are not reported. 221 of 225 metabolites had a coefficient of variation below 0.1 across 14 cell lines with biological replicates. About 9% of values were missing and were imputed with MICE.

**Culture and harvest.** Each line was grown in its own (vendor-recommended) medium, which the authors list in their Supplementary Table 1. Cells were seeded at 2 million per T75 flask and extracted after 2 days, before reaching 90% confluence. Polar metabolites were extracted with cold 80% methanol.

**Main findings relevant to this project.**
- Lineage is a strong driver of metabolite levels: 148 of 225 metabolites had lineage r² > 0.1 (e.g., 1-methylnicotinamide is high in most kidney and pleural lines).
- Some of the clearest associations are **genetic or epigenetic rather than transcriptional**: IDH1/2 hotspot mutations and 2-hydroxyglutarate; KEAP1 mutation and GSH/GSSG and NADP⁺; ME2 deletion and malate; ASNS promoter hypermethylation and asparaginase sensitivity.
- Several enzyme-level links were shown with expression or methylation: SLC25A20 and long-chain acylcarnitines, PYCR1 and proline, GPT2 and alanine, IDO1/IDO2/TDO2 and kynurenine.
- Common oncogenes and tumour suppressors (EGFR, KRAS, NRAS, TP53, PTEN, TSC1/2) showed weak or no metabolite associations.

**Implications for this project.**
1. Values are **relative** (normalised peak areas). They can be compared across cell lines for one metabolite, not between metabolites.
2. Lineage and medium must be treated as baseline predictors, not nuisances to be ignored.
3. Metabolites whose variation is driven by mutations or deletions (2HG, malate) are expected to be poorly explained by expression of their network genes. This is useful for interpretation, not a flaw.
4. Imputation was done by the data generators on the full dataset. This cannot be undone, and it is a small, unavoidable source of information sharing between samples.
5. There is a reporting discrepancy downstream: Cavicchioli et al. (2022) describe the panel as 136 polar and 89 non-polar metabolites, whereas Li et al. (2019) state 124 polar and 101 lipid species. The project should classify metabolites from the assay method annotation in the original supplementary data, not from secondary papers.

---

## 3. Studies that directly predicted metabolite levels from expression in CCLE

### 3.1 Li, Barbour, Zhu & Wong (2022), *FASEB J* — "poor predictor"

- **Data.** 454 solid, adherent cell lines (184 blood cancer lines excluded, and only adherent lines kept to reduce culture-method bias); 225 metabolites; RSEM expression filtered to 13,403 genes; 18 media types (RPMI 1640 49.7%, DMEM 15.0%, EMEM 11.7%).
- **Medium and clinical covariates alone.** In stepwise regression with 55 medium components and six clinical traits, only 7 metabolites reached R² > 0.5, and 206 of 225 (91.6%) had R² < 0.3. Strong multicollinearity was noted.
- **Model.** Single random split (318 train / 136 test). Genes were pre-ranked by univariable R² and the top 10% kept, then LASSO (glmnet, 10-fold CV for λ). The paper does not clearly state whether the univariable pre-filter was restricted to the training set.
- **Results.** Median test R² was about 0 (range −0.52 to 0.69); 115 of 225 metabolites (51%) had test R² ≤ 0. Training R² had a median of 0.58, indicating heavy overfitting. Only a handful of metabolites exceeded R² 0.4 on the test set; 1-methylnicotinamide (MNA), thymidine and phosphocreatine are named. NNMT was the top predictor of MNA.
- **Metabolic vs all genes.** Restricting to 554 SMPDB or 1,701 Reactome metabolic genes gave median test R² near 0, generally lower than all-gene models. Pathway-level GSVA scores computed from genes and from metabolites correlated weakly (−0.18 to 0.25).
- **Bulk tissue.** In a breast cancer tissue dataset (n = 102), NNMT–MNA showed no significant relationship, which the authors partly attribute to tumour purity.
- **Critical reading.** A single 70/30 split gives high-variance estimates for each metabolite. Thresholds for "accurate" are used inconsistently (R² > 0.4 in test; R² > 0.5 in both sets). LASSO is unstable with thousands of correlated predictors, and a ridge or elastic-net model might have done somewhat better. Pathway membership (SMPDB/Reactome) is much coarser than reaction-level neighbourhoods. Even so, the central message — that held-out R² is low for most metabolites — is plausible and should be treated as the prior.

### 3.2 Cavicchioli et al. (2022), *Int J Mol Sci* — "can be predicted"

- **Data.** CCLE expression (1,019 lines) and metabolomics (928 lines), overlap of about 890–898 lines; NCI-60 (57 lines) for cross-dataset tests.
- **Method.** Correlation networks between each metabolite and transcripts (Pearson, p < 0.001), pruned with the Data Processing Inequality and bootstrapped (the `corto` package). Predictions are network-weighted enrichment scores, i.e., relative values rather than calibrated abundances.
- **Validation.** Repeated 50/50 splits within CCLE. Across 157 metabolites (those with at least one significant direct transcript), the mean predicted–measured correlation was 0.424, and all 157 had positive mean correlation. Degree-preserving shuffled networks gave a mean of 0.00. Best within-CCLE examples include MNA, arginine, fructose-1-phosphate, valine and pyroglutamic acid. Transfer between CCLE and NCI-60 worked for a subset of shared metabolites (e.g., methionine, tyrosine, leucine, phenylalanine, uracil); allantoin transferred negatively in both directions.
- **Enzyme links noted.** NNMT–MNA, IDO1–kynurenine, DPYD–N-carbamoyl-β-alanine, and SLC6A8 (creatine transporter) with creatine and phosphocreatine. Networks also contained many non-metabolic genes (e.g., MAX, ATM, PTEN, BRAF).
- **Critical reading.** Lineage and medium were **not** controlled. A correlation of 0.42 corresponds to roughly 18% of variance if predictions were perfectly calibrated, so much of this signal may be lineage structure that any expression-based predictor would capture. The shuffled-network null shows the network is better than nothing; it does not show it is better than lineage, or better than a random set of expressed genes of the same size. There are internal inconsistencies in the number of splits reported (1,000 vs 2,000). TCGA "predictions" were tested only for survival association, not against measured metabolites.

### 3.3 Schildcrout et al. (2024 preprint "Recon8D"; published 2026 as "Metab8D", *Commun Biol*)

- **Data and model.** CCLE, eight omics classes (CNV, DNA methylation, histone PTMs, transcriptomics, splicing, miRNA, proteomics, phosphoproteomics) predicting the 225-metabolite panel. Random forest was chosen after benchmarking against ridge, LASSO, XGBoost and SVM. 5-fold CV, Pearson correlation, Bonferroni significance.
- **Findings (preprint).** Transcriptomics gave the most accurate predictions of all omics classes. Restricting to fewer than 2,000 Recon3D metabolic genes still predicted over 99% of metabolites significantly, with a slightly lower median correlation than all genes. **Models built on 2,000 random genes performed similarly**, so the metabolic restriction was not what drove accuracy. Among top transcriptomic features, the authors report that about 17% were enzymes in the metabolite's own reactions and most of the rest were in first-degree neighbouring reactions (calculated over 95 overlapping metabolites; this figure should be checked in the final article). Known pairs were recovered (NNMT–MNA, SLC6A6–taurine), alongside signalling genes such as YAP1.
- **Confounding.** Controls stratified by growth type, one lineage contrast and tumour origin. Culture medium was not modelled as a covariate.
- **Critical reading.** "Significantly predicted" with n ≈ 900 requires only a modest correlation, so 96–99% significance does not mean accurate prediction. The random-gene result is the single most important prior finding for this project: it suggests that a large share of the predictive signal reflects general cell state (lineage, growth rate, differentiation), not local metabolic control.

### 3.4 How these three results fit together

The studies appear to disagree ("poor predictor" vs "can be predicted" vs "tightly associated"), but most of the disagreement is explained by **what was measured and what was controlled**:

| Study | Cell lines | Model | Validation | Lineage/medium controlled? | Metric | Headline |
|---|---|---|---|---|---|---|
| Li 2022 | 454 solid adherent | Univariable filter + LASSO | One 70/30 split | Partly (blood/suspension lines removed; medium and lineage also tried in models, but not as a formal baseline) | Test R² | Median test R² ≈ 0 |
| Cavicchioli 2022 | ~890 all | Correlation network (`corto`) | Repeated 50/50 splits; NCI-60 | No | Pearson r | Mean r ≈ 0.42 |
| Recon8D 2024 | ~900 (per supp. table) | Random forest | 5-fold CV | Stratified checks; medium not modelled | Pearson r, significance | Transcriptome best; random genes ≈ metabolic genes |

Three reconciling points:
1. **Correlation vs calibrated R².** A held-out Pearson r of 0.4 and a held-out R² near 0.1–0.2 can describe the same predictions. R² can also be negative when predictions are miscalibrated even if they are correlated.
2. **Confounding.** Without a lineage-plus-medium baseline, correlations partly measure how well expression encodes tissue of origin.
3. **Feature specificity.** None of the three studies asked, per metabolite, whether the genes that make or consume that metabolite carry more information than random genes of the same size after accounting for lineage and medium. Recon8D came closest, but compared large gene pools (2,000 genes), not metabolite-specific neighbourhoods.

---

## 4. Evidence on network proximity and gene–metabolite covariation

### 4.1 Zelezniak, Sheridan & Patil (2014), *PLoS Comput Biol*
In yeast, the authors found that kinetic constraints at single reactions described mRNA–metabolite relationships poorly, while constraints from network connectivity (the production and consumption reactions surrounding a metabolite) were key determinants. This is the conceptual basis for using producing and consuming genes as predictors. **Caveat:** it was shown in one organism across dynamic perturbations, not across genetically diverse human cell lines in different media.

### 4.2 Benedetti et al. (2023), *Nature Metabolism*
- A meta-analysis of 988 tumour and normal specimens across 11 cancer types (15 datasets), using **Human1** to define network distance between gene and metabolite (distance 1 = shared reaction; 2 = consecutive reactions; "proximal" = 1 or 2).
- Only about 10.7% of gene–metabolite pairs had a computable network distance.
- Of 22,619 significant gene–metabolite associations, only 565 (about 2.5%) were proximal. Proximal pairs were enriched relative to non-significant pairs (odds ratio 1.42), but the authors describe proximity as a weak determinant of the overall association landscape.
- A second class of associations involved hub metabolites (e.g., quinolinate, NAD⁺-related metabolites) correlating with immune-cell genes, i.e., cell composition rather than pathway control.
- In about 900 CCLE lines, IDO1 expression associated with kynurenine, while AFMID (the next enzyme) did not.
- **Relevance.** This is the closest precedent for using Human1 network distance in gene–metabolite analysis. It predicts modest but real enrichment of network-local signal. Its tissue setting (bulk tumours, mixed cell types) differs from cell lines, where composition effects are much smaller — which arguably makes cell lines a cleaner test of the local-control hypothesis.

### 4.3 Ortmayr, Dubuis & Zampieri (2019), *Nature Communications*
Profiling of 54 adherent NCI-60 lines with transcriptome and proteome integration uncovered widespread associations between transcriptional regulators and metabolic pathways. This supports the alternative hypothesis that metabolite levels are shaped by **global regulatory programmes**, not only by local enzyme abundance. (Reviewed at abstract level.)

### 4.4 Siddiqui et al. (2018), *BMC Bioinformatics* (IntLIM)
A linear-model framework (metabolite ~ gene + phenotype + gene × phenotype) for phenotype-specific gene–metabolite associations, applied to NCI-60. Relevant mainly as precedent for simple, interpretable linear models in this setting.

---

## 5. Metabolic gene expression structure in tumours (Peng et al., 2018)

Using 9,125 TCGA samples across 33 cancer types, Peng et al. defined tumour subtypes from mRNA expression of seven major metabolic processes. Subtypes correlated with survival (upregulated carbohydrate, nucleotide and vitamin/cofactor metabolism with worse prognosis; upregulated lipid metabolism with better prognosis), with somatic drivers, and with master regulators; SNAI1 or RUNX1 knockdown altered metabolic activity and drug sensitivity.

**Relevance.** This study did not measure metabolite abundance at scale, so it does not test expression → metabolite prediction. Its importance here is indirect but significant: **metabolic gene expression is strongly structured by tumour type and regulatory programme.** Any set of metabolic genes will therefore carry lineage and cell-state information. This is exactly why the project needs (a) a lineage-plus-medium baseline and (b) a null of random *metabolic* genes, not only random genes from the whole genome.

---

## 6. Culture medium as a confounder

- **Vande Voorde et al. (2019), *Science Advances*.** Showed that commonly used media (e.g., DMEM) impose metabolic artefacts on cancer cells and that a physiological medium (Plasmax) changes cellular metabolism.
- **Joly, Chew & Graham (2021), *PLoS Comput Biol*.** In CCLE/DepMap analyses of metabolic pathway dependencies, culture medium was confirmed as a major confounding factor.
- **Li et al. (2022).** Medium and clinical traits alone explained little variance for most metabolites (91.6% with R² < 0.3), but medium is collinear with lineage (e.g., haematopoietic lines are mostly in RPMI).
- **Cherkaoui et al. (2022), *Mol Syst Biol*.** 180 cancer cell lines grown in a single medium (phenol-red-free RPMI 1640, 2 mM glutamine, 10% FCS, at least two weeks of adaptation), profiled by untargeted flow-injection MS (1,809 ions putatively annotated; no chromatographic separation, so isomers are not resolved; annotation broadly at confidence level 4). The authors argue that sequencing alone does not capture pathway activity or carbon flux. This dataset removes medium confounding but has much lower identification confidence and a different panel; it is **not** recommended for the initial project (see `research_plan.md`, Section 12).

**Implication.** Intracellular levels of nutrients supplied by the medium (e.g., essential amino acids, glucose, glutamine) are expected to reflect medium composition and uptake more than enzyme expression. These targets should be interpreted separately.

---

## 7. A methodological warning from gene-signature research (Venet et al., 2011)

Venet, Dumont and Detours showed that 28 of 47 published breast cancer outcome signatures (60%) were not significantly better than random signatures of the same size, and that over 90% of random signatures with more than 100 genes were significantly associated with outcome. The reason was that most random gene sets capture a dominant biological axis (proliferation). Adjusting for that axis removed most of the association.

**Relevance.** The same trap applies here. A gene set that "predicts" a metabolite is meaningful only if it beats size-matched random gene sets processed through the identical pipeline. Recon8D's random-gene result (Section 3.3) suggests this trap is real for CCLE metabolomics.

---

## 8. Tools: Human1/Human-GEM and COBRApy

### 8.1 Human1 (Robinson et al., 2020, *Science Signaling*)
- A consensus human genome-scale metabolic model integrating HMR2, iHsa and Recon3D. At publication it had 13,417 reactions, 10,138 metabolites (4,164 unique, the rest being compartment-specific copies) and 3,625 genes, and achieved 100% stoichiometric consistency in Memote tests (compared with 19.8% for Recon3D).
- Distributed through a version-controlled GitHub repository (SysBioChalmers/Human-GEM) with SBML, YAML and MATLAB files, plus annotation tables (`reactions.tsv`, `metabolites.tsv`, `genes.tsv`) holding external identifiers. Genes are identified by Ensembl gene IDs.

### 8.2 The model has changed substantially since 2020
- Release notes document repeated gene and GPR curation: removal of 558 genes attached to generic (de)phosphorylation reactions (v1.11), 149 non-metabolic genes (v1.14), 21 pseudogenes (v1.15), HMDB IDs updated to version 4 format (v1.13), and compartment abbreviation changes (v1.11).
- **Luo et al. (2026), *PNAS*** describe "Human2" (Human-GEM v2.0.0): an LLM-assisted review of 26,246 gene–reaction pairs led to 1,135 GPR updates, removal of 203 genes unrelated to metabolism and 717 revised lipid GPRs. The resulting model has 2,848 genes, 12,931 reactions and 8,461 metabolites. The authors explicitly note that human GEMs contain incorrect GPR associations and reactions added for gap-filling or simulation purposes.
- **Implication.** Gene sets derived from Human-GEM depend on the version. The project must pin one release (by tag and commit hash) and should not mix annotations across versions. A cheap robustness check is to compare target gene sets between v1.19 and v2.0.

### 8.3 COBRApy (Ebrahim et al., 2013, *BMC Systems Biology*)
COBRApy reads SBML models and exposes reactions, metabolites, genes and gene–protein–reaction (GPR) rules as Python objects (e.g., each metabolite's reactions, each reaction's genes and its rule string). The Human-GEM README lists COBRApy among tools that can use its SBML file; the repository's own helper functions are MATLAB-only.

**Assessment.** COBRApy is appropriate for this project, but only as a **parser and graph-query tool**. No flux balance analysis, solver or simulation is needed. Pandas reading of the annotation TSV files is needed alongside it to map metabolites to HMDB/KEGG/ChEBI identifiers. The scientific risks lie in the model content (GPR quality, hub metabolites, compartments, identifier mapping), not in the software. These are discussed in `research_plan.md`, Sections 3 and 10.

---

## 9. Synthesis

### 9.1 What is reasonably established
1. Gene expression contains statistically detectable information about many CCLE metabolite levels (Cavicchioli 2022; Recon8D).
2. For most metabolites, held-out predictive accuracy is low once evaluated as calibrated R² (Li 2022).
3. A few metabolites have strong, reproducible single-enzyme relationships: MNA–NNMT (Li 2022, Cavicchioli 2022, Recon8D), kynurenine–IDO1/TDO2 (Li 2019, Cavicchioli 2022, Benedetti 2023), creatine/phosphocreatine–SLC6A8 (Cavicchioli 2022), N-carbamoyl-β-alanine–DPYD (Cavicchioli 2022), taurine–SLC6A6 (Recon8D).
4. Lineage is a major source of metabolite variation (Li 2019), and culture medium is a known confounder (Joly 2021; Vande Voorde 2019).
5. Several strong metabolite phenotypes are driven by mutation, deletion or methylation (2HG, malate, asparagine), not necessarily captured by expression of network genes (Li 2019).
6. Network proximity is enriched among gene–metabolite associations in tumours but explains only a small fraction of them (Benedetti 2023).

### 9.2 Where findings conflict or remain unclear
1. **Whether prediction is "good" or "poor"** — largely a matter of metric (correlation vs R²) and confounder control (Section 3.4).
2. **Whether metabolic genes are special.** Li 2022 found metabolic gene subsets did no better (often worse) than all genes; Recon8D found metabolic genes about equal to random genes; Recon8D also reported that top features were enriched near the metabolite in the network. These are compatible only if a small number of local genes carry specific signal for a few metabolites, while most predictive power is diffuse cell-state information.
3. **Whether the network signal survives a lineage-plus-medium baseline.** Not directly tested by any of the three CCLE prediction studies.

### 9.3 What remains worth investigating
As far as this review found, no study has performed the following specific, controlled comparison:

> For each of a set of reliably identified metabolites, do the genes encoding enzymes and transporters of reactions that produce or consume that metabolite (from a curated genome-scale model) improve held-out prediction **beyond lineage and culture medium**, and do they do so **more than size-matched random gene sets** — including random *metabolic* gene sets — evaluated with the same cross-validation folds?

This is a narrow question. It is **not** novel to show that gene–metabolite associations exist, nor to predict CCLE metabolites from expression. The contribution of this project would be a careful, well-controlled test of a specific mechanistic hypothesis on a public dataset, with an honest answer that may well be negative or mixed. That is appropriate for a portfolio project: the value lies in sound design and interpretation, not in a claim of discovery.

---

## 10. References

1. Li H, Ning S, Ghandi M, et al. The landscape of cancer cell line metabolism. *Nat Med.* 2019;25(5):850–860. doi:10.1038/s41591-019-0404-8
2. Li H, Barbour JA, Zhu X, Wong JWH. Gene expression is a poor predictor of steady-state metabolite abundance in cancer cells. *FASEB J.* 2022;36(5):e22296. doi:10.1096/fj.202101921RR (preprint: bioRxiv doi:10.1101/2021.12.19.473333)
3. Cavicchioli MV, Santorsola M, Balboni N, Mercatelli D, Giorgi FM. Prediction of metabolic profiles from transcriptomics data in human cancer cell lines. *Int J Mol Sci.* 2022;23(7):3867. doi:10.3390/ijms23073867
4. Peng X, Chen Z, Farshidfar F, et al. Molecular characterization and clinical relevance of metabolic expression subtypes in human cancers. *Cell Rep.* 2018;23(1):255–269.e4. doi:10.1016/j.celrep.2018.03.077
5. Schildcrout R, Smith K, Bhowmick R, et al. Recon8D: a metabolic regulome network from oct-omics and machine learning. *bioRxiv.* 2024. doi:10.1101/2024.08.17.608400. Published as: Metab8D: a metabolic regulome network from multiomics and machine learning. *Commun Biol.* 2026. doi:10.1038/s42003-026-10409-3
6. Benedetti E, Liu EM, Tang C, et al. A multimodal atlas of tumour metabolism reveals the architecture of gene–metabolite covariation. *Nat Metab.* 2023;5(6):1029–1044. doi:10.1038/s42255-023-00817-8
7. Zelezniak A, Sheridan S, Patil KR. Contribution of network connectivity in determining the relationship between gene expression and metabolite concentration changes. *PLoS Comput Biol.* 2014;10(4):e1003572. doi:10.1371/journal.pcbi.1003572
8. Ortmayr K, Dubuis S, Zampieri M. Metabolic profiling of cancer cells reveals genome-wide crosstalk between transcriptional regulators and metabolism. *Nat Commun.* 2019;10:1841. doi:10.1038/s41467-019-09695-9
9. Siddiqui JK, Baskin E, Liu M, et al. IntLIM: integration using linear models of metabolomics and gene expression data. *BMC Bioinformatics.* 2018. doi:10.1186/s12859-018-2085-6
10. Vande Voorde J, Ackermann T, Pfetzer N, et al. Improving the metabolic fidelity of cancer models with a physiological cell culture medium. *Sci Adv.* 2019;5(1):eaau7314. doi:10.1126/sciadv.aau7314
11. Joly JH, Chew BTL, Graham NA. The landscape of metabolic pathway dependencies in cancer cell lines. *PLoS Comput Biol.* 2021;17(4):e1008942. doi:10.1371/journal.pcbi.1008942
12. Cherkaoui S, Durot S, Bradley J, et al. A functional analysis of 180 cancer cell lines reveals conserved intrinsic metabolic programs. *Mol Syst Biol.* 2022;18:e11033. doi:10.15252/msb.202211033
13. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol.* 2011;7(10):e1002240. doi:10.1371/journal.pcbi.1002240
14. Robinson JL, Kocabaş P, Wang H, et al. An atlas of human metabolism. *Sci Signal.* 2020;13(624):eaaz1482. doi:10.1126/scisignal.aaz1482
15. Luo J, Wang H, Moyer D, et al. Reconstruction of human metabolic models with large language models. *Proc Natl Acad Sci USA.* 2026;123(15):e2516511123.
16. Ebrahim A, Lerman JA, Palsson BO, Hyduke DR. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. *BMC Syst Biol.* 2013;7:74. doi:10.1186/1752-0509-7-74
17. Ghandi M, Huang FW, Jané-Valbuena J, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature.* 2019;569:503–508. doi:10.1038/s41586-019-1186-3
18. Brunk E, Sahoo S, Zielinski DC, et al. Recon3D enables a three-dimensional view of gene variation in human metabolism. *Nat Biotechnol.* 2018;36:272–281. doi:10.1038/nbt.4072
19. SysBioChalmers. Human-GEM repository and release notes. https://github.com/SysBioChalmers/Human-GEM (accessed 16 November 2026).
