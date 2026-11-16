"""Phase 3 compute: main analysis and the four pre-specified sensitivity analyses.

    .venv/bin/python -m src.phase3.run main     # primary set: M0, M1, M1-enz, M2, N-met, N-all, P  (74 targets + MNA diagnostic)
    .venv/bin/python -m src.phase3.run sens2    # primary set: N-met size-matched to the enzymatic-only set (67 targets)
    .venv/bin/python -m src.phase3.run sens3    # primary set: variance-matched N-met
    .venv/bin/python -m src.phase3.run sens1    # all cancer lines incl. haematopoietic (896): M0, M1, N-met
    .venv/bin/python -m src.phase3.run sens4    # primary minus Cellosaurus-problematic lines (709): M0, M1, N-met

Each target's arrays are written to results/phase3/<analysis>/<feature>.npz once; existing files are skipped,
so an interrupted run resumes. Every analysis uses B = 1,000 (D45) and the frozen gene sets.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd

from ..phase0.sources import REPO
from ..phase2 import folds as fo
from ..phase2 import pipeline as pl
from ..phase2.data import load_inputs
from ..phase2.run import FOLDS as PRIMARY_FOLDS, run_target

OUT = REPO / "results" / "phase3"
B = 1000
FOLD_SEEDS = {"all_lines": 2000, "no_problematic": 3000}       # D46; primary folds are the frozen Phase 2 folds

ANALYSES = {
    #  name    population        kwargs for run_target
    "main":  ("primary",        dict(with_m2=True, with_nall=True)),
    "sens2": ("primary",        dict(with_m2=False, with_nall=False, with_enz_null=True)),
    "sens3": ("primary",        dict(with_m2=False, with_nall=False, with_varmatched=True)),
    "sens1": ("all_lines",      dict(with_m2=False, with_nall=False)),
    "sens4": ("no_problematic", dict(with_m2=False, with_nall=False)),
}


def _safe(name: str) -> str:
    return name.replace("/", "_").replace(" ", "_")


def folds_for(population: str, cells: pd.DataFrame) -> pd.DataFrame:
    if population == "primary":
        f = pd.read_csv(PRIMARY_FOLDS)
    else:
        path = OUT / f"folds_{population}.csv"
        if path.exists():
            f = pd.read_csv(path)
        else:
            f = fo.make_folds(cells, seed0=FOLD_SEEDS[population])
            f.to_csv(path, index=False, lineterminator="\n")
    fo.validate_folds(f, cells)
    return f


def run_analysis(analysis: str, only: list[str] | None = None) -> None:
    population, kw = ANALYSES[analysis]
    inp = load_inputs(population)
    folds = folds_for(population, inp.cells)
    ev = pl.Evaluator(inp, fo.fold_indices(folds, inp.cells))
    var_dec = pl.variance_deciles(inp.X, inp.idx(inp.met_pool)) if kw.get("with_varmatched") else None
    out = OUT / analysis
    out.mkdir(parents=True, exist_ok=True)
    feats = inp.targets + (inp.diagnostics if analysis == "main" else [])
    if only:
        feats = [f for f in feats if f in only]
    meta = {"analysis": analysis, "population": population, "n_lines": int(len(inp.cells)), "B": B,
            "folds": "results/phase2/folds.csv" if population == "primary" else f"results/phase3/folds_{population}.csv",
            "n_features": len(feats), "kwargs": {k: v for k, v in kw.items()}, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    t_all = time.time()
    for i, name in enumerate(feats):
        path = out / f"{_safe(name)}.npz"
        if path.exists():
            continue
        t = time.time()
        res = run_target(name, inp, ev, B=B, var_dec=var_dec, **kw)
        arrays = {k: v for k, v in res.items() if isinstance(v, np.ndarray)}
        scalars = {k: v for k, v in res.items() if not isinstance(v, np.ndarray)}
        np.savez_compressed(path, **arrays, _scalars=json.dumps(scalars))
        print(f"[{analysis}] {i + 1:3d}/{len(feats)} {name:32s} k={res['k']:3d} {time.time() - t:6.1f}s", flush=True)
    print(f"[{analysis}] done in {(time.time() - t_all) / 60:.1f} min", flush=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("analysis", choices=list(ANALYSES))
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args(argv)
    run_analysis(a.analysis, a.only)
    return 0


if __name__ == "__main__":
    sys.exit(main())
