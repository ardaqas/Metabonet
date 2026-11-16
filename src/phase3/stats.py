"""Phase 3 statistics: per-target metrics, global tests and hypothesis classification exactly as fixed in
expectations.md (§1 definitions, §3 H1-H6, D36/D37). Reads results/phase3/<analysis>/*.npz.

    .venv/bin/python -m src.phase3.stats
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu, wilcoxon

from ..phase0.sources import PROCESSED, REPO
from ..phase2 import folds as fo
from ..phase2 import pipeline as pl
from ..phase2.data import load_inputs
from ..phase2.run import FOLDS as PRIMARY_FOLDS

OUT = REPO / "results" / "phase3"
P1 = PROCESSED / "phase1"
MIN_EFFECT = 0.02          # expectations §1
FDR_Q = 0.10               # expectations H2
ALPHA = 0.05
CAPTURE_MIN_M2 = 0.02      # research_plan §9


def _safe(name: str) -> str:
    return name.replace("/", "_").replace(" ", "_")


def bh(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q-values (NaN-aware)."""
    q = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    pv = p[ok].to_numpy()
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    qv = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(qv, 1.0)
    q[ok] = out
    return q


def emp_p(null_draws: np.ndarray, observed: float) -> float:
    """(1 + #{draws >= observed}) / (B + 1); draws are per-draw means over repeats."""
    return (1 + int((null_draws >= observed).sum())) / (len(null_draws) + 1)


def load_analysis(analysis: str) -> tuple[pd.DataFrame, dict]:
    """Per-target table plus per-draw null means (for figures)."""
    sel = pd.read_csv(P1 / "target_selection.csv").set_index("feature")
    clusters = pd.read_csv(P1 / "target_clusters.csv").set_index("feature")
    d = OUT / analysis
    rows, nulls = [], {}
    for path in sorted(d.glob("*.npz")):
        z = np.load(path)
        sc = json.loads(str(z["_scalars"])) if "_scalars" in z.files else {}
        name = sc["feature"]
        r2m0 = z["r2_m0"]
        d1 = z["r2_m1"] - r2m0
        rec = {"feature": name, "k": int(sc["k"]), "k_enz": int(sc.get("k_enz", 0)),
               "class": sel.at[name, "class"] if name in sel.index else "",
               "positive_control": bool(sel.at[name, "positive_control"]) if name in sel.index else False,
               "role": "diagnostic" if bool(sel.at[name, "pipeline_diagnostic"]) else "target",
               "cluster_id": clusters["cluster_id"].get(name, ""),
               "r2_m0": r2m0.mean(), "rho_m0": z["rho_m0"].mean(),
               "dr2_m1": d1.mean(), "dr2_m1_min_rep": d1.min(), "dr2_m1_max_rep": d1.max(), "rho_m1": z["rho_m1"].mean(),
               "dr2_perm": (z["r2_m1_perm"] - z["r2_m0_perm"]).mean()}
        if "r2_m2" in z.files:
            d2 = z["r2_m2"] - r2m0
            rec.update(dr2_m2=d2.mean(), rho_m2=z["rho_m2"].mean())
            rec["capture_ratio"] = rec["dr2_m1"] / rec["dr2_m2"] if rec["dr2_m2"] > CAPTURE_MIN_M2 else np.nan
        if "r2_m1_enz" in z.files:
            rec["dr2_m1_enz"] = (z["r2_m1_enz"] - r2m0).mean()
        for key, lab in (("d_nmet", "nmet"), ("d_nall", "nall"), ("d_nmet_enz", "nmet_enz"), ("d_nmet_vm", "nmet_vm")):
            if key in z.files:
                dm = z[key].mean(axis=1)                       # per-draw mean over repeats
                nulls.setdefault(lab, {})[name] = dm
                obs = rec["dr2_m1_enz"] if lab == "nmet_enz" else rec["dr2_m1"]
                rec[f"{lab}_median"] = float(np.median(dm))
                rec[f"{lab}_p95"] = float(np.percentile(dm, 95))
                rec[f"p_{lab}"] = emp_p(dm, obs)
                rec[f"B_{lab}"] = len(dm)
        rows.append(rec)
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df, nulls
    tg = df["role"] == "target"
    for lab in ("nmet", "nall", "nmet_enz", "nmet_vm"):
        if f"p_{lab}" in df:
            df[f"q_{lab}"] = np.nan
            df.loc[tg, f"q_{lab}"] = bh(df.loc[tg, f"p_{lab}"])
            obs = df["dr2_m1_enz"] if lab == "nmet_enz" else df["dr2_m1"]
            df[f"network_specific_{lab}"] = tg & (df[f"q_{lab}"] < FDR_Q) & (obs >= MIN_EFFECT)
            df[f"above_p95_{lab}"] = tg & (obs > df[f"{lab}_p95"])
    return df.sort_values("feature").reset_index(drop=True), nulls


# ------------------------------------------------------------------------------------------------- global tests
def _wilcoxon_gt0(x: pd.Series) -> float:
    x = x.dropna()
    x = x[x != 0]
    return float(wilcoxon(x, alternative="greater").pvalue) if len(x) >= 5 else np.nan


def _binom(n_above: int, n: int, p0: float = 0.05) -> float:
    return float(binomtest(n_above, n, p0, alternative="greater").pvalue) if n else np.nan


def h2_class(global_p: float, n_specific: int) -> str:
    """expectations H2: supported / partially supported / not supported."""
    if global_p < ALPHA and n_specific >= 3:
        return "supported"
    if n_specific == 0 and not (global_p < ALPHA):
        return "not supported"
    if 1 <= n_specific <= 2 and not (global_p < ALPHA):
        return "partially supported"
    # cases the pre-specification does not name explicitly (D49): reported as partial with the reason
    if n_specific >= 3:
        return "partially supported (>=3 network-specific targets but global p >= 0.05)"
    return "partially supported (global p < 0.05 but < 3 network-specific targets)"


def null_test_block(t: pd.DataFrame, lab: str, obs_col: str = "dr2_m1") -> dict:
    """H2/H3-style block against one null family."""
    diff = t[obs_col] - t[f"{lab}_median"]
    n_spec = int(t[f"network_specific_{lab}"].sum())
    n_above = int(t[f"above_p95_{lab}"].sum())
    p_glob = _wilcoxon_gt0(diff)
    return {"n_targets": int(len(t)), "wilcoxon_p_vs_median_null": p_glob, "median_excess_over_null": float(diff.median()),
            "n_above_own_p95": n_above, "expected_above_p95": 0.05 * len(t), "binomial_p": _binom(n_above, len(t)),
            "n_network_specific": n_spec, "network_specific_targets": sorted(t.loc[t[f"network_specific_{lab}"], "feature"]),
            "classification": h2_class(p_glob, n_spec)}


def cluster_block(t: pd.DataFrame, lab: str = "nmet") -> dict:
    """D37 secondary reading: one observation per Jaccard cluster."""
    g = t.assign(excess=t["dr2_m1"] - t[f"{lab}_median"], above=t["dr2_m1"] - t[f"{lab}_p95"]).groupby("cluster_id")
    per = g.agg(n=("feature", "size"), excess=("excess", "median"), above=("above", "median")).reset_index()
    return {"n_clusters": int(len(per)), "n_multi_target_clusters": int((per["n"] > 1).sum()),
            "wilcoxon_p": _wilcoxon_gt0(per["excess"]), "median_excess": float(per["excess"].median()),
            "n_clusters_above_p95": int((per["above"] > 0).sum()), "binomial_p": _binom(int((per["above"] > 0).sum()), len(per)),
            "per_cluster": per.to_dict("records")}


def global_tests(df: pd.DataFrame) -> dict:
    t = df[df["role"] == "target"].copy()
    out = {"n_targets": int(len(t))}
    # S1 on the main run
    out["S1_permuted_control"] = {"median": float(t["dr2_perm"].median()), "max": float(t["dr2_perm"].max()),
                                  "n_above_0.02": int((t["dr2_perm"] > 0.02).sum()),
                                  "pass": bool(abs(t["dr2_perm"].median()) <= 0.005 and (t["dr2_perm"] <= 0.02).all())}
    # H1
    p1 = _wilcoxon_gt0(t["dr2_m1"])
    n_ge = int((t["dr2_m1"] >= MIN_EFFECT).sum())
    out["H1"] = {"wilcoxon_p": p1, "median_dr2_m1": float(t["dr2_m1"].median()), "n_dr2_ge_0.02": n_ge,
                 "frac_dr2_ge_0.02": n_ge / len(t), "n_dr2_positive": int((t["dr2_m1"] > 0).sum()),
                 "supported": bool(p1 < ALPHA and n_ge / len(t) >= 0.25)}
    # H2 (primary: N-met), secondary clustered reading, H3 (N-all)
    out["H2"] = null_test_block(t, "nmet")
    out["H2_clustered_secondary"] = cluster_block(t, "nmet")
    out["H3"] = null_test_block(t, "nall")
    # H4
    if "dr2_m2" in t:
        d42 = t["dr2_m2"] - t["dr2_m1"]
        cap = t["capture_ratio"].dropna()
        out["H4"] = {"wilcoxon_p_m2_gt_m1": _wilcoxon_gt0(d42), "median_dr2_m2": float(t["dr2_m2"].median()),
                     "n_m2_ge_m1": int((d42 >= 0).sum()), "n_m1_beats_m2_by_0.02": int((-d42 > 0.02).sum()),
                     "m1_beats_m2_targets": sorted(t.loc[-d42 > 0.02, "feature"]),
                     "n_capture_defined": int(len(cap)), "median_capture_ratio": float(cap.median()) if len(cap) else np.nan,
                     "n_capture_ge_0.7": int((cap >= 0.7).sum()), "capture_ge_0.7_targets": sorted(t.loc[t["capture_ratio"] >= 0.7, "feature"]),
                     "supported": bool(_wilcoxon_gt0(d42) < ALPHA)}
    # H5
    pcs = t[t["positive_control"]]
    out["H5"] = {"controls": [{"feature": r.feature, "dr2_m1": r.dr2_m1, "p_nmet": r.p_nmet, "q_nmet": r.q_nmet,
                               "passes": bool(r.dr2_m1 >= MIN_EFFECT and r.p_nmet < ALPHA)} for r in pcs.itertuples()],
                 "n_pass": int(sum(r.dr2_m1 >= MIN_EFFECT and r.p_nmet < ALPHA for r in pcs.itertuples())), "n_controls": int(len(pcs))}
    # H6 (exploratory)
    ex = t["dr2_m1"] - t["nmet_median"]
    cls = {c: t[t["class"] == c] for c in ("synthesised", "medium-supplied", "mixed")}
    mw = mannwhitneyu(ex[t["class"] == "medium-supplied"], ex[t["class"] == "synthesised"], alternative="two-sided") \
        if (t["class"] == "medium-supplied").sum() >= 3 else None
    out["H6"] = {"by_class": {c: {"n": int(len(g)), "median_r2_m0": float(g["r2_m0"].median()), "median_dr2_m1": float(g["dr2_m1"].median()),
                                  "median_excess_over_nmet": float((g["dr2_m1"] - g["nmet_median"]).median()),
                                  "n_network_specific": int(g["network_specific_nmet"].sum())} for c, g in cls.items() if len(g)},
                 "mannwhitney_p_medium_vs_synth": float(mw.pvalue) if mw else np.nan}
    # diagnostic
    dg = df[df["role"] == "diagnostic"]
    out["pipeline_diagnostic"] = [{"feature": r.feature, "k": r.k, "dr2_m1": r.dr2_m1, "p_nmet": r.p_nmet, "nmet_median": r.nmet_median,
                                   "note": "not a target; outside selection rules (D34)"} for r in dg.itertuples()]
    return out


def unadjusted_controls(df: pd.DataFrame) -> pd.DataFrame:
    """H5 diagnostic: R^2 of the gene set alone (no lineage/medium stage), for controls + diagnostic (D48)."""
    inp = load_inputs("primary")
    ev = pl.Evaluator(inp, fo.fold_indices(pd.read_csv(PRIMARY_FOLDS), inp.cells))
    feats = list(df.loc[df["positive_control"] | (df["role"] == "diagnostic"), "feature"])
    rows = []
    for name in feats:
        y = inp.Y[name].to_numpy(dtype=float)
        tf = ev.target_folds(y)
        # unadjusted: baseline = training mean only
        base_te = [[np.full(len(te), y[tr].mean()) for tr, te in reps] for reps in ev.folds]
        res_tr = [[y[tr] - y[tr].mean() for tr, _ in reps] for reps in ev.folds]
        sst = ((y - y.mean()) ** 2).sum()
        r2_mean = []
        for reps, b in zip(ev.folds, base_te):
            pred = np.empty(len(y))
            for (tr, te), p in zip(reps, b):
                pred[te] = p
            r2_mean.append(1 - ((y - pred) ** 2).sum() / sst)
        tf0 = pl.TargetFolds(base_te, res_tr, np.array(r2_mean), np.zeros(len(ev.folds)))
        g = inp.idx(inp.sets_primary[name])
        m1u = ev.evaluate(y, tf0, g)
        rows.append({"feature": name, "r2_genes_only_unadjusted": float((m1u.r2 - tf0.r2_m0).mean()),
                     "r2_m0": float(tf.r2_m0.mean()), "dr2_m1_adjusted": float(df.set_index("feature").at[name, "dr2_m1"])})
    return pd.DataFrame(rows)


def sensitivity_summary(main: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """H2 classification under each sensitivity analysis (expectations §6; D36)."""
    rows, blocks = [], {}
    mt = main[main["role"] == "target"]
    blocks["main"] = {"analysis": "main (primary set, N-met)", **{k: v for k, v in null_test_block(mt, "nmet").items() if k != "network_specific_targets"}}
    for ana, lab, obs, label in (("sens1", "nmet", "dr2_m1", "1: all cancer lines incl. haematopoietic (896)"),
                                 ("sens2", "nmet_enz", "dr2_m1_enz", "2: enzymatic-only gene sets (targets with k_enz >= 3)"),
                                 ("sens3", "nmet_vm", "dr2_m1", "3: variance-matched random metabolic genes"),
                                 ("sens4", "nmet", "dr2_m1", "4: primary set without 26 problematic lines (709)")):
        df, _ = load_analysis(ana)
        if len(df) == 0:
            continue
        t = df[df["role"] == "target"]
        if ana == "sens2":
            t = t[t["k_enz"] >= 3].copy()
            t[f"q_{lab}"] = bh(t[f"p_{lab}"])      # BH within the 67 analysed targets
            t[f"network_specific_{lab}"] = (t[f"q_{lab}"] < FDR_Q) & (t[obs] >= MIN_EFFECT)
        blk = null_test_block(t, lab, obs)
        blocks[ana] = {"analysis": label, **{k: v for k, v in blk.items() if k != "network_specific_targets"},
                       "network_specific_targets": blk["network_specific_targets"]}
        for r in t.itertuples():
            rows.append({"analysis": ana, "feature": r.feature, "dr2": getattr(r, obs), f"p_null": getattr(r, f"p_{lab}"),
                         "q_null": getattr(r, f"q_{lab}"), "network_specific": bool(getattr(r, f"network_specific_{lab}"))})
    classes = {k: v["classification"] for k, v in blocks.items()}
    base = classes["main"].split(" (")[0]
    blocks["robustness"] = {"H2_main": classes["main"], "H2_by_analysis": classes,
                            "robust": all(c.split(" (")[0] == base for c in classes.values()),
                            "n_analyses_available": len(classes) - 1}
    return pd.DataFrame(rows), blocks


def main() -> None:
    df, nulls = load_analysis("main")
    df.to_csv(OUT / "per_target_main.csv", index=False, float_format="%.6g", lineterminator="\n")
    np.savez_compressed(OUT / "null_draw_means_main.npz", **{f"{lab}__{_safe(n)}": v for lab, d in nulls.items() for n, v in d.items()})
    g = global_tests(df)
    uc = unadjusted_controls(df)
    uc.to_csv(OUT / "h5_controls_unadjusted.csv", index=False, float_format="%.6g", lineterminator="\n")
    g["H5"]["unadjusted"] = uc.to_dict("records")
    srows, sblocks = sensitivity_summary(df)
    srows.to_csv(OUT / "per_target_sensitivity.csv", index=False, float_format="%.6g", lineterminator="\n")
    g["sensitivity"] = sblocks
    (OUT / "global_tests.json").write_text(json.dumps(g, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
    print(json.dumps({k: g[k] for k in ("H1", "H2", "H3") if k in g}, indent=1, default=str)[:3000])


if __name__ == "__main__":
    main()
