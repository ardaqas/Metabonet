# MetaboNet

Does the expression of a metabolite's Human-GEM first-neighbour genes predict its abundance in CCLE cell lines beyond lineage and culture medium, and better than size-matched random gene sets? The specification is in `docs/research_plan.md`, `docs/expectations.md` and `docs/literature_review.md`.

**Status:** Phase 0 (data audit, `docs/data_audit_report.md`), Phase 1 (network, mapping, frozen 74-target list, `docs/phase1_report.md`) Phase 2 (CV pipeline and controls, `docs/phase2_report.md`) and Phase 3 (main analysis and sensitivity analyses, `docs/phase3_report.md`) are complete. Phase 4 (case studies and write-up) has not started.

```
data/raw/        downloaded files (not committed; provenance in data/raw/MANIFEST.tsv)
data/processed/  frozen Phase 0 tables + FREEZE.tsv hashes
src/phase0/      download, Cellosaurus, loaders, medium normalisation, crosswalk, audits
src/phase1/      Human-GEM parsing, gene/metabolite mapping, gene sets, target selection
src/phase2/      folds, closed-form ridge, two-stage CV pipeline, controls (run_target = Phase 3 entry point)
src/phase3/      main + sensitivity runs (per-target npz), statistics (H1–H6), figures, report tables
results/         folds and control outputs (phase2/), later metrics and null distributions
notebooks/       00_data_audit.ipynb (executed)
tests/           unit tests
docs/            plan, expectations, literature review, decision log, data sources, audit report
```

Reproduce: `make env && make data && make test && make audit && make phase1 && make phase2 && make phase3`. `data/raw/depmap/Model.csv` and `README.txt` must be downloaded manually, and `RELEASE.txt` written by hand; see `docs/data_sources.md`.
