PY := .venv/bin/python
NB := notebooks/00_data_audit.ipynb

.PHONY: env data verify freeze audit test all phase1 phase2 phase3

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

all: data test audit

phase1:         ## Phase 1: Human-GEM parsing, gene + metabolite mapping, gene sets, outcome-blind target selection
	$(PY) -m src.phase0.download --verify
	$(PY) -m src.phase1.run

phase2:         ## Phase 2: frozen folds + pipeline controls (S1 permuted target, S2 synthetic recovery, S5 timing)
	$(PY) -m src.phase0.download --verify
	$(PY) -m src.phase2.run all

phase3:         ## Phase 3: main analysis + 4 sensitivity analyses (resumable; ~3 h total), statistics, figures, report tables
	$(PY) -m src.phase0.download --verify
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run main
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens2
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens3
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens1
	OMP_NUM_THREADS=2 $(PY) -m src.phase3.run sens4
	$(PY) -m src.phase3.manifest
	$(PY) -m src.phase3.stats
	$(PY) -m src.phase3.figures
	$(PY) -m src.phase3.report_tables
