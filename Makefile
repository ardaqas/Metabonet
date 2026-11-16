PY := .venv/bin/python
NB := notebooks/00_data_audit.ipynb

.PHONY: env data verify freeze audit test all phase1 phase2 phase3 phase4 verify-results check-docs

env:            ## create the pinned virtual environment
	uv venv --python 3.10 .venv && uv pip install --python $(PY) -r requirements.txt

data:           ## download raw data (skips existing files), query Cellosaurus, write data/raw/MANIFEST.tsv
	$(PY) -m src.phase0.download
	$(PY) -m src.phase0.cellosaurus

verify:         ## re-hash raw files and compare with MANIFEST (no network)
	$(PY) -m src.phase0.download --verify

freeze: verify  ## build and freeze data/processed/ tables without the notebook
	$(PY) -m src.phase0.run

audit: verify   ## execute the audit notebook in place (it rebuilds + freezes data/processed/)
	$(PY) -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=1200 $(NB)

test:
	$(PY) -m pytest -q tests

all: data test audit phase1 phase2 phase3 phase4   ## whole project from raw data (manual step: Model.csv + README.txt + RELEASE.txt in data/raw/depmap/)

phase1:         ## Phase 1: Human-GEM parsing, gene + metabolite mapping, gene sets, outcome-blind target selection
	$(PY) -m src.phase0.download --verify
	$(PY) -m src.phase1.run

phase2:         ## Phase 2: frozen folds + pipeline controls (S1 permuted target, S2 synthetic recovery, S5 timing)
	$(PY) -m src.phase0.download --verify
	$(PY) -m src.phase2.run all
	$(PY) -m src.phase2.run s2x

phase3:         ## Phase 3: main + 4 sensitivity analyses (resumable; ~3 h), statistics, post hoc other-target null, figures, tables
	$(PY) -m src.phase0.download --verify
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run main
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens2
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens3
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens1
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens4
	$(PY) -m src.phase3.manifest
	$(PY) -m src.phase3.stats
	$(PY) -m src.phase3.other_null
	$(PY) -m src.phase3.figures
	$(PY) -m src.phase3.report_tables

phase4:         ## Phase 4: case studies (gene roles, coefficients, out-of-fold predictions) and case-study figures
	$(PY) -m src.phase4.casestudies
	$(PY) -m src.phase4.figures

verify-results: ## check the Phase 3 per-target arrays against results/phase3/RESULTS_MANIFEST.tsv
	$(PY) -m src.phase3.manifest --verify

check-docs:     ## check relative links, absolute paths and file sizes, and the numbers quoted in README.md
	$(PY) scripts/check_docs.py all
