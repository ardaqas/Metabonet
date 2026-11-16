# MetaboNet

Does the expression of a metabolite's Human-GEM first-neighbour genes predict its abundance in CCLE cell lines beyond lineage and culture medium, and better than size-matched random gene sets? The specification is in `docs/research_plan.md`, `docs/expectations.md` and `docs/literature_review.md`.

**Status:** Phase 0 (data audit) is complete; see `docs/data_audit_report.md`. Phase 1 has not started.

```
data/raw/        downloaded files (not committed; provenance in data/raw/MANIFEST.tsv)
data/processed/  frozen Phase 0 tables + FREEZE.tsv hashes
src/phase0/      download, Cellosaurus, loaders, medium normalisation, crosswalk, audits
notebooks/       00_data_audit.ipynb (executed)
tests/           unit tests
docs/            plan, expectations, literature review, decision log, data sources, audit report
```

Reproduce: `make env && make data && make test && make audit`. `data/raw/depmap/Model.csv` and `README.txt` must be downloaded manually, and `RELEASE.txt` written by hand; see `docs/data_sources.md`.
