"""Phase 2: pipeline controls (expectations.md §2).

    .venv/bin/python -m src.phase2.run folds   # generate + save folds (once)
    .venv/bin/python -m src.phase2.run s1      # permuted-target control, all 74 targets + diagnostic
    .venv/bin/python -m src.phase2.run s2      # synthetic-target recovery
    .venv/bin/python -m src.phase2.run s5      # time one complete target (D38/D45)
    .venv/bin/python -m src.phase2.run all

`run_target` is the single entry point that Phase 3 will reuse unchanged.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd

from ..phase0.sources import REPO
from . import folds as fo
from . import pipeline as pl
from .data import load_inputs, strata

OUT = REPO / "results" / "phase2"
FOLDS = OUT / "folds.csv"
B_DEFAULT = 1000
S2_N_SIMS, S2_N_SIGNAL_GENES, S2_B, S2_SIGNAL_VAR = 20, 5, 200, 0.15   # D44


def get_evaluator():
    inp = load_inputs()
    if not FOLDS.exists():
        raise FileNotFoundError("run `python -m src.phase2.run folds` first")
    f = pd.read_csv(FOLDS)
    fo.validate_folds(f, inp.cells)
    return inp, pl.Evaluator(inp, fo.fold_indices(f, inp.cells))


def cmd_folds() -> None:
    inp = load_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    if FOLDS.exists():
        print(f"{FOLDS} exists; not regenerated (folds are frozen)")
        return
    f = fo.make_folds(inp.cells)
    f.to_csv(FOLDS, index=False, lineterminator="\n")
    print(f"wrote {FOLDS} ({f.repeat.nunique()} repeats x {f.fold.nunique()} folds)")


def run_target(name: str, inp, ev: pl.Evaluator, B: int = B_DEFAULT, with_m2: bool = True,
               with_enzymatic: bool = True, with_nall: bool = True, with_enz_null: bool = False,
               with_varmatched: bool = False, var_dec: np.ndarray | None = None) -> dict:
    """Full evaluation of one feature under the frozen folds. Returns arrays only (no interpretation)."""
    y = inp.Y[name].to_numpy(dtype=float)
    tf = ev.target_folds(y)
    g = inp.idx(inp.sets_primary[name])
    k = len(g)
    out = {"feature": name, "k": k, "r2_m0": tf.r2_m0, "rho_m0": tf.rho_m0}
    t = time.time()
    m1 = ev.evaluate(y, tf, g)
    out.update(r2_m1=m1.r2, rho_m1=m1.rho, alpha_m1=m1.alpha_idx, t_m1=time.time() - t)
    if with_enzymatic:
        ge = inp.idx(inp.sets_enzymatic[name])
        if len(ge) >= 1:
            me = ev.evaluate(y, tf, ge)
            out.update(k_enz=len(ge), r2_m1_enz=me.r2, rho_m1_enz=me.rho)
    if with_m2:
        t = time.time()
        m2 = ev.evaluate(y, tf, np.arange(inp.X.shape[1]))
        out.update(r2_m2=m2.r2, rho_m2=m2.rho, alpha_m2=m2.alpha_idx, t_m2=time.time() - t)
    if B > 0:
        met_pool = np.setdiff1d(inp.idx(inp.met_pool), g)
        all_pool = np.setdiff1d(np.arange(inp.X.shape[1]), g)
        t = time.time()
        d_met, dr_met = ev.null_draws(y, tf, met_pool, k, B, np.random.default_rng(pl.target_seed(name, "nmet")))
        out.update(d_nmet=d_met, draws_nmet=dr_met, t_nmet=time.time() - t)
        if with_nall:
            t = time.time()
            d_all, dr_all = ev.null_draws(y, tf, all_pool, k, B, np.random.default_rng(pl.target_seed(name, "nall")))
            out.update(d_nall=d_all, draws_nall=dr_all, t_nall=time.time() - t)
        if with_enz_null and out.get("k_enz", 0) >= 3:        # sensitivity 2 (D36): size-matched to k_enz
            ge = inp.idx(inp.sets_enzymatic[name])
            pool_e = np.setdiff1d(inp.idx(inp.met_pool), ge)
            d_e, dr_e = ev.null_draws(y, tf, pool_e, len(ge), B, np.random.default_rng(pl.target_seed(name, "nmet_enz")))
            out.update(d_nmet_enz=d_e, draws_nmet_enz=dr_e)
        if with_varmatched:                                      # sensitivity 3 (D47): variance-decile matched
            assert var_dec is not None
            dr_v = pl.matched_draws(met_pool, var_dec, g, B, np.random.default_rng(pl.target_seed(name, "nmet_vm")))
            d_v, _ = ev.null_draws(y, tf, met_pool, k, B, None, draws=dr_v)
            out.update(d_nmet_vm=d_v, draws_nmet_vm=dr_v)
    # permuted-target control (D43): one stratified permutation, identical folds and procedure
    t = time.time()
    yp = pl.permute_within_strata(y, strata(inp.cells), np.random.default_rng(pl.target_seed(name, "perm")))
    tfp = ev.target_folds(yp)
    mp = ev.evaluate(yp, tfp, g)
    out.update(r2_m0_perm=tfp.r2_m0, r2_m1_perm=mp.r2, t_perm=time.time() - t)
    return out


def cmd_s1() -> None:
    inp, ev = get_evaluator()
    st = strata(inp.cells)
    rows = []
    for name in inp.targets + inp.diagnostics:
        y = inp.Y[name].to_numpy(dtype=float)
        g = inp.idx(inp.sets_primary[name])
        yp = pl.permute_within_strata(y, st, np.random.default_rng(pl.target_seed(name, "perm")))
        tfp = ev.target_folds(yp)
        mp = ev.evaluate(yp, tfp, g)
        d = mp.r2 - tfp.r2_m0
        rows.append({"feature": name, "role": "diagnostic" if name in inp.diagnostics else "target", "k": len(g),
                     "r2_m0_perm": tfp.r2_m0.mean(), "r2_m1_perm": mp.r2.mean(), "delta_r2_perm": d.mean(),
                     "delta_r2_perm_min": d.min(), "delta_r2_perm_max": d.max(), "delta_r2_perm_sd": d.std(ddof=1)})
        print(f"  {name:32s} k={len(g):3d}  dR2_P={d.mean():+.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s1_permuted_control.csv", index=False, float_format="%.6g", lineterminator="\n")
    tg = df[df["role"] == "target"]
    res = {"n_targets": int(len(tg)), "median_delta_r2_perm": float(tg["delta_r2_perm"].median()),
           "max_delta_r2_perm": float(tg["delta_r2_perm"].max()), "min_delta_r2_perm": float(tg["delta_r2_perm"].min()),
           "n_targets_above_0.02": int((tg["delta_r2_perm"] > 0.02).sum()),
           "pass_median_within_0.005": bool(abs(tg["delta_r2_perm"].median()) <= 0.005),
           "pass_no_target_above_0.02": bool((tg["delta_r2_perm"] <= 0.02).all())}
    res["S1_pass"] = res["pass_median_within_0.005"] and res["pass_no_target_above_0.02"]
    (OUT / "s1_summary.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))


def cmd_s2() -> None:
    """Synthetic-target recovery (D44): y = lineage + medium effects + 5-gene linear signal + noise."""
    inp, ev = get_evaluator()
    rng0 = np.random.default_rng(20261116)
    cands = sorted([f for f in inp.targets if len(inp.sets_primary[f]) >= S2_N_SIGNAL_GENES],
                   key=lambda f: (len(inp.sets_primary[f]), f))
    step = max(1, len(cands) // S2_N_SIMS)
    sims = cands[::step][:S2_N_SIMS]
    lin = pd.factorize(inp.cells["lineage_collapsed"])[0]
    med = pd.factorize(inp.cells["medium_category"])[0]
    Xz = (inp.X.astype(float) - inp.X.mean(0)) / np.where(inp.X.std(0) == 0, 1, inp.X.std(0))
    rows = []
    for i, name in enumerate(sims):
        rng = np.random.default_rng(rng0.integers(2**31))
        g = inp.idx(inp.sets_primary[name])
        sig = rng.choice(g, size=S2_N_SIGNAL_GENES, replace=False)
        beta = rng.normal(size=S2_N_SIGNAL_GENES)
        s = Xz[:, sig] @ beta
        s *= np.sqrt(S2_SIGNAL_VAR) / s.std()
        le = rng.normal(size=lin.max() + 1)[lin]; le *= np.sqrt(0.30) / le.std()
        mde = rng.normal(size=med.max() + 1)[med]; mde *= np.sqrt(0.05) / mde.std()
        eps = rng.normal(size=len(s)); eps *= np.sqrt(0.50) / eps.std()
        y = le + mde + s + eps
        true_delta = s.var() / y.var()
        tf = ev.target_folds(y)
        m1 = ev.evaluate(y, tf, g)
        d1 = m1.delta(tf)
        pool = np.setdiff1d(inp.idx(inp.met_pool), g)
        dn, _ = ev.null_draws(y, tf, pool, len(g), S2_B, np.random.default_rng(rng.integers(2**31)))
        p = (1 + (dn.mean(axis=1) >= d1.mean()).sum()) / (S2_B + 1)
        rows.append({"sim": i, "gene_set_of": name, "k": len(g), "signal_genes": ";".join(inp.genes[sig]),
                     "true_delta_r2": true_delta, "r2_m0": tf.r2_m0.mean(), "delta_r2_m1": d1.mean(),
                     "delta_r2_m1_sd_repeats": d1.std(ddof=1), "recovery_ratio": d1.mean() / true_delta,
                     "nmet_median_delta": float(np.median(dn.mean(axis=1))), "p_vs_nmet": p, "beats_nmet": p < 0.05})
        print(f"  sim {i:2d} {name:28s} k={len(g):2d} true={true_delta:.3f} M1={d1.mean():.3f} p={p:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s2_synthetic.csv", index=False, float_format="%.6g", lineterminator="\n")
    res = {"n_sims": len(df), "B_nmet": S2_B, "signal_var_target": S2_SIGNAL_VAR,
           "median_true_delta": float(df["true_delta_r2"].median()), "median_delta_m1": float(df["delta_r2_m1"].median()),
           "median_recovery_ratio": float(df["recovery_ratio"].median()),
           "median_abs_error": float((df["delta_r2_m1"] - df["true_delta_r2"]).abs().median()),
           "frac_beats_nmet": float(df["beats_nmet"].mean()),
           "pass_beats_nmet_ge_95pct": bool(df["beats_nmet"].mean() >= 0.95),
           "pass_close_abs_error_le_0.05": bool((df["delta_r2_m1"] - df["true_delta_r2"]).abs().median() <= 0.05)}
    res["S2_pass"] = res["pass_beats_nmet_ge_95pct"] and res["pass_close_abs_error_le_0.05"]
    (OUT / "s2_summary.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))


def cmd_s2_extended(n_sims: int = 200) -> None:
    """S2 rerun with `n_sims` simulations (D53). Same design as D44 (signal share 0.15, 5 signal genes, B = 200),
    fresh seeds, cycling through every target with >= 5 genes; written to separate files so the original
    20-simulation result (s2_synthetic.csv) is untouched."""
    inp, ev = get_evaluator()
    rng0 = np.random.default_rng(20261117)
    cands = sorted([f for f in inp.targets if len(inp.sets_primary[f]) >= S2_N_SIGNAL_GENES],
                   key=lambda f: (len(inp.sets_primary[f]), f))
    lin = pd.factorize(inp.cells["lineage_collapsed"])[0]
    med = pd.factorize(inp.cells["medium_category"])[0]
    Xz = (inp.X.astype(float) - inp.X.mean(0)) / np.where(inp.X.std(0) == 0, 1, inp.X.std(0))
    rows = []
    for i in range(n_sims):
        name = cands[i % len(cands)]
        rng = np.random.default_rng(rng0.integers(2**31))
        g = inp.idx(inp.sets_primary[name])
        sig = rng.choice(g, size=S2_N_SIGNAL_GENES, replace=False)
        beta = rng.normal(size=S2_N_SIGNAL_GENES)
        sgl = Xz[:, sig] @ beta; sgl *= np.sqrt(S2_SIGNAL_VAR) / sgl.std()
        le = rng.normal(size=lin.max() + 1)[lin]; le *= np.sqrt(0.30) / le.std()
        mde = rng.normal(size=med.max() + 1)[med]; mde *= np.sqrt(0.05) / mde.std()
        eps = rng.normal(size=len(sgl)); eps *= np.sqrt(0.50) / eps.std()
        y = le + mde + sgl + eps
        true_delta = sgl.var() / y.var()
        tf = ev.target_folds(y)
        d1 = ev.evaluate(y, tf, g).delta(tf)
        pool = np.setdiff1d(inp.idx(inp.met_pool), g)
        dn, _ = ev.null_draws(y, tf, pool, len(g), S2_B, np.random.default_rng(rng.integers(2**31)))
        p = (1 + (dn.mean(axis=1) >= d1.mean()).sum()) / (S2_B + 1)
        rows.append({"sim": i, "gene_set_of": name, "k": len(g), "true_delta_r2": true_delta, "r2_m0": tf.r2_m0.mean(),
                     "delta_r2_m1": d1.mean(), "recovery_ratio": d1.mean() / true_delta,
                     "nmet_median_delta": float(np.median(dn.mean(axis=1))), "p_vs_nmet": p, "beats_nmet": p < 0.05})
        if (i + 1) % 20 == 0:
            print(f"  {i + 1}/{n_sims} sims; beat N-met so far {np.mean([r['beats_nmet'] for r in rows]):.3f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s2_200_synthetic.csv", index=False, float_format="%.6g", lineterminator="\n")
    err = (df["delta_r2_m1"] - df["true_delta_r2"]).abs()
    res = {"n_sims": len(df), "n_distinct_gene_sets": int(df["gene_set_of"].nunique()), "B_nmet": S2_B,
           "median_true_delta": float(df["true_delta_r2"].median()), "median_delta_m1": float(df["delta_r2_m1"].median()),
           "median_recovery_ratio": float(df["recovery_ratio"].median()),
           "recovery_ratio_p05_p95": [float(df["recovery_ratio"].quantile(.05)), float(df["recovery_ratio"].quantile(.95))],
           "median_abs_error": float(err.median()), "frac_beats_nmet": float(df["beats_nmet"].mean()),
           "n_fail_beats_nmet": int((~df["beats_nmet"]).sum()),
           "pass_beats_nmet_ge_95pct": bool(df["beats_nmet"].mean() >= 0.95),
           "pass_close_abs_error_le_0.05": bool(err.median() <= 0.05)}
    res["S2_200_pass"] = res["pass_beats_nmet_ge_95pct"] and res["pass_close_abs_error_le_0.05"]
    (OUT / "s2_200_summary.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))


def cmd_s5(B: int = B_DEFAULT) -> None:
    """Time one complete target (D45): the target with the median primary gene-set size."""
    inp, ev = get_evaluator()
    order = sorted(inp.targets, key=lambda f: (len(inp.sets_primary[f]), f))
    name = order[len(order) // 2]
    t0 = time.time()
    res = run_target(name, inp, ev, B=B)
    wall = time.time() - t0
    ks = np.array([len(inp.sets_primary[f]) for f in inp.targets])
    k0 = res["k"]
    # null cost scales roughly with k (SVD of n x k); project per-target cost with a k ratio on the null part
    null_t = res["t_nmet"] + res["t_nall"]
    fixed_t = wall - null_t
    projected = float(np.sum(fixed_t + null_t * np.maximum(ks, 3) / max(k0, 3)))
    summary = {"timing_target": name, "k": int(k0), "B": B, "wall_seconds": wall,
               "t_m1": res["t_m1"], "t_m2": res["t_m2"], "t_nmet": res["t_nmet"], "t_nall": res["t_nall"], "t_perm": res["t_perm"],
               "projected_total_hours_74_targets": projected / 3600,
               "contingency_threshold_hours": 8.0, "reduce_B_to_500": projected / 3600 > 8.0,
               "note": "result arrays stored in timing_target_*.npz; not interpreted in Phase 2"}
    np.savez_compressed(OUT / f"timing_target_{name.replace('/', '_')}.npz",
                        **{k: v for k, v in res.items() if isinstance(v, np.ndarray)})
    (OUT / "s5_timing.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["folds", "s1", "s2", "s2x", "s5", "all"])
    ap.add_argument("--B", type=int, default=B_DEFAULT)
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.cmd in ("folds", "all"):
        cmd_folds()
    if a.cmd in ("s1", "all"):
        cmd_s1()
    if a.cmd in ("s2", "all"):
        cmd_s2()
    if a.cmd == "s2x":
        cmd_s2_extended()
    if a.cmd in ("s5", "all"):
        cmd_s5(a.B)
    return 0


if __name__ == "__main__":
    sys.exit(main())
