# MetaboNet

Does the expression of a metabolite's Human-GEM first-neighbour genes predict its abundance in CCLE cell lines beyond lineage and culture medium, and better than size-matched random gene sets? The specification is in `docs/research_plan.md`, `docs/expectations.md` and `docs/literature_review.md`.

**Status:** all four phases are complete. Start with **[`docs/report.md`](docs/report.md)** (the final write-up), then [`docs/case_studies.md`](docs/case_studies.md). Phase reports: [data audit](docs/data_audit_report.md), [network and targets](docs/phase1_report.md), [pipeline and controls](docs/phase2_report.md), [main results](docs/phase3_report.md). Every decision: [`docs/decision_log.md`](docs/decision_log.md); the hypotheses fixed in advance: [`docs/expectations.md`](docs/expectations.md).

**Result in one line:** for 74 pre-specified polar metabolites in 735 solid-tumour cell lines, the genes in each metabolite's own reactions predict its held-out abundance beyond lineage and medium and beyond size-matched random and other-pathway gene sets (37 of 74 targets individually; robust to four sensitivity analyses), with strong heterogeneity and the limitations listed in the report.

```
data/raw/        downloaded files (not committed; provenance in data/raw/MANIFEST.tsv)
data/processed/  frozen Phase 0 tables + FREEZE.tsv hashes
src/phase0/      download, Cellosaurus, loaders, medium normalisation, crosswalk, audits
src/phase1/      Human-GEM parsing, gene/metabolite mapping, gene sets, target selection
src/phase2/      folds, closed-form ridge, two-stage CV pipeline, controls (run_target = Phase 3 entry point)
src/phase3/      main + sensitivity runs (per-target npz), statistics (H1–H6), post hoc other-target null, figures, tables
src/phase4/      case studies (gene roles, coefficients, out-of-fold predictions) and their figures
results/         phase2/ controls and folds, phase3/ per-target results, statistics and figures, phase4/ case studies
notebooks/       00_data_audit.ipynb (executed)
tests/           unit tests
docs/            plan, expectations, literature review, decision log, data sources, audit report
```

Reproduce: `make env && make data && make all` (about 5 h; thread count is fixed for bit-identical results). `make verify-results` checks the stored Phase 3 arrays against `results/phase3/RESULTS_MANIFEST.tsv`. `data/raw/depmap/Model.csv` and `README.txt` must be downloaded manually, and `RELEASE.txt` written by hand; see `docs/data_sources.md`.
