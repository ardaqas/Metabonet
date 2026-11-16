"""POST HOC exploratory analysis (D54): other-target null.

Added after seeing the Phase 3 results; does not change any pre-specified conclusion.

For each target m: fit the frozen primary expressed gene sets of every *other* target j (j != m, MNA diagnostic
excluded) that shares no gene with G_m and has 0.5k <= |G_j| <= 2k (k = |G_m|) to m's own metabolite, with
run_target() on the frozen folds, and compare m's own Delta-R^2(M1) against that null.

    .venv/bin/python -m src.phase3.other_null
Writes results/phase3/other_target_null.csv, other_target_null_sets.csv, other_target_null_summary.json.
"""
from __future__ import annotations

import dataclasses
import json
import time

import numpy as np
import pandas as pd

from ..phase0.sources import PROCESSED, REPO
from ..phase2 import folds as fo
from ..phase2 import pipeline as pl
from ..phase2.data import load_inputs
from ..phase2.run import FOLDS, run_target

OUT = REPO / "results" / "phase3"
MIN_ELIGIBLE = 10                         # D54 (4)
AMINO_ACIDS = ["alanine", "arginine", "asparagine", "aspartate", "glutamate", "glutamine", "glycine", "histidine",
               "isoleucine", "leucine", "lysine", "methionine", "phenylalanine", "proline", "serine", "threonine",
               "tryptophan", "tyrosine", "valine"]


def group_of(feature: str, cluster_id: str, cluster_size: int) -> str:
    if cluster_size > 1 and cluster_id == "C01":
        return "nucleosides"
    if cluster_size > 1 and cluster_id == "C02":
        return "carnitines"
    if feature in AMINO_ACIDS:
        return "amino acids"
    return "others"


def eligible_sets(sets: dict[str, list[str]], targets: list[str], m: str) -> list[str]:
    """Other targets whose frozen set is disjoint from G_m and within 0.5k..2k (D54 rule 1)."""
    gm, k = set(sets[m]), len(sets[m])
    return [j for j in targets if j != m and not (set(sets[j]) & gm) and 0.5 * k <= len(sets[j]) <= 2 * k]


def fit_other(inp, ev, m: str, j: str) -> float:
    """Delta-R^2 of target m's metabolite using the gene set of target j, via run_target() (D54 rule 2)."""
    swapped = dataclasses.replace(inp, sets_primary={**inp.sets_primary, m: inp.sets_primary[j]})
    res = run_target(m, swapped, ev, B=0, with_m2=False, with_enzymatic=False)
    return float((res["r2_m1"] - res["r2_m0"]).mean())


def main() -> None:
    inp = load_inputs("primary")
    ev = pl.Evaluator(inp, fo.fold_indices(pd.read_csv(FOLDS), inp.cells))
    main_df = pd.read_csv(OUT / "per_target_main.csv").set_index("feature")
    cl = pd.read_csv(PROCESSED / "phase1" / "target_clusters.csv").set_index("feature")
    targets = list(inp.targets)

    # sanity: the swapped-set path reproduces Phase 3 for a target's own set (identical code path)
    for m in ("creatine", "2-deoxycytidine"):
        own = fit_other(inp, ev, m, m)
        z = np.load(OUT / "main" / f"{m}.npz")                      # exact arrays (the CSV keeps 6 digits)
        exact = float((z["r2_m1"] - z["r2_m0"]).mean())
        assert abs(own - exact) < 1e-12, (m, own, exact)

    t0 = time.time()
    rows, detail = [], []
    for i, m in enumerate(targets):
        elig = eligible_sets(inp.sets_primary, targets, m)
        d = np.array([fit_other(inp, ev, m, j) for j in elig])
        for j, v in zip(elig, d):
            detail.append({"feature": m, "other_set_of": j, "k_other": len(inp.sets_primary[j]), "dr2": v})
        obs = float(main_df.at[m, "dr2_m1"])
        n = len(d)
        rec = {"feature": m, "class": main_df.at[m, "class"], "group": group_of(m, cl.at[m, "cluster_id"], int(cl.at[m, "cluster_size"])),
               "cluster_id": cl.at[m, "cluster_id"], "k": len(inp.sets_primary[m]), "n_eligible": n,
               "too_few_to_judge": n < MIN_ELIGIBLE, "dr2_m1": obs,
               "network_specific_nmet": bool(main_df.at[m, "network_specific_nmet"]), "nmet_median": float(main_df.at[m, "nmet_median"])}
        if n:
            rec.update(other_median=float(np.median(d)), other_p95=float(np.percentile(d, 95)), other_max=float(d.max()),
                       percentile_rank=100.0 * float((d < obs).sum()) / n, emp_p=(1 + int((d >= obs).sum())) / (n + 1),
                       exceeds_p95=bool(obs > np.percentile(d, 95)))
        rows.append(rec)
        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{len(targets)} targets, {time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    det = pd.DataFrame(detail)
    df.to_csv(OUT / "other_target_null.csv", index=False, float_format="%.6g", lineterminator="\n")
    det.to_csv(OUT / "other_target_null_sets.csv", index=False, float_format="%.6g", lineterminator="\n")

    ns = df[df["network_specific_nmet"]]
    judged = ns[~ns["too_few_to_judge"]]
    win = judged[judged["exceeds_p95"]]
    frac = len(win) / len(judged) if len(judged) else float("nan")
    reading = ("specific neighbourhood" if frac > 0.5 else ("coherent metabolic modules" if frac < 0.5 else "mixed"))

    def block(g: pd.DataFrame) -> dict:
        j = g[~g["too_few_to_judge"]]
        return {"n_network_specific": int(len(g)), "n_judged": int(len(j)), "n_exceed_p95": int(j["exceeds_p95"].sum()),
                "n_too_few": int(g["too_few_to_judge"].sum()),
                "median_percentile_rank": float(j["percentile_rank"].median()) if len(j) else None}
    by_group = {grp: block(ns[ns["group"] == grp]) for grp in ["nucleosides", "carnitines", "amino acids", "others"]}
    allt = df[~df["too_few_to_judge"]]
    summary = {
        "label": "POST HOC, added after seeing results; does not change any pre-specified conclusion (D54)",
        "n_targets": int(len(df)), "n_network_specific": int(len(ns)),
        "eligible_sets_per_target": {"min": int(df["n_eligible"].min()), "median": float(df["n_eligible"].median()), "max": int(df["n_eligible"].max())},
        "n_targets_too_few": int(df["too_few_to_judge"].sum()),
        "network_specific_too_few": sorted(ns.loc[ns["too_few_to_judge"], "feature"]),
        "network_specific_judged": int(len(judged)), "network_specific_exceed_p95": int(len(win)),
        "fraction_exceed_p95": frac, "reading_rule_outcome": reading,
        "network_specific_not_exceeding": sorted(judged.loc[~judged["exceeds_p95"], "feature"]),
        "by_group": by_group,
        "all_targets_judged_exceed_p95": int(allt["exceeds_p95"].sum()), "all_targets_judged": int(len(allt)),
        "median_of_other_medians": float(df["other_median"].median()), "median_of_nmet_medians": float(df["nmet_median"].median()),
        "median_other_p95": float(df["other_p95"].median()), "median_percentile_rank_network_specific": float(judged["percentile_rank"].median()),
        "n_non_network_specific_exceed_p95": int(df[~df["network_specific_nmet"] & ~df["too_few_to_judge"]]["exceeds_p95"].sum()),
    }
    # resolution and margin diagnostics (reported alongside the pre-stated criterion)
    jj = judged.assign(floor=1.0 / (judged["n_eligible"] + 1))
    at_floor = jj[(jj["emp_p"] >= 0.05) & np.isclose(jj["emp_p"], jj["floor"])]
    summary["resolution"] = {
        "n_emp_p_lt_0.05": int((jj["emp_p"] < 0.05).sum()),
        "n_emp_p_ge_0.05": int((jj["emp_p"] >= 0.05).sum()),
        "emp_p_ge_0.05_at_resolution_floor_(exceed_every_comparator)": sorted(at_floor["feature"]),
        "emp_p_ge_0.05_not_at_floor": sorted(set(jj.loc[jj["emp_p"] >= 0.05, "feature"]) - set(at_floor["feature"])),
        "n_rank_100": int((jj["percentile_rank"] == 100).sum()),
        "below_best_comparator": [{"feature": r.feature, "n_eligible": int(r.n_eligible), "dr2_m1": r.dr2_m1, "best_comparator": r.other_max,
                                   "emp_p": r.emp_p} for r in jj[jj["dr2_m1"] <= jj["other_max"]].itertuples()],
        "min_n_eligible_judged_network_specific": int(jj["n_eligible"].min()),
        "median_margin_over_best_comparator": float((jj["dr2_m1"] - jj["other_max"]).median()),
    }
    summary["pooled_other_target_fits"] = {
        "n_fits": int(len(det)), "median": float(det["dr2"].median()), "p95": float(det["dr2"].quantile(.95)),
        "p99": float(det["dr2"].quantile(.99)), "max": float(det["dr2"].max()), "frac_ge_0.02": float((det["dr2"] >= 0.02).mean()),
        "top5": det.nlargest(5, "dr2").round(4).to_dict("records")}
    summary["judged_not_network_specific_exceeding_p95"] = sorted(
        df[~df["network_specific_nmet"] & ~df["too_few_to_judge"] & df["exceeds_p95"]]["feature"])
    summary["too_few_network_specific_descriptive"] = [
        {"feature": r.feature, "n_eligible": int(r.n_eligible), "dr2_m1": r.dr2_m1, "best_comparator": r.other_max}
        for r in df[df["network_specific_nmet"] & df["too_few_to_judge"]].itertuples()]
    (OUT / "other_target_null_summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
