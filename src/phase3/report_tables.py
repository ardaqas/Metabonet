"""Render Phase 3 markdown tables into docs/phase3_report.md from docs/phase3_report.template.md.

    .venv/bin/python -m src.phase3.report_tables
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from ..phase0.sources import REPO

OUT = REPO / "results" / "phase3"


def _md(df: pd.DataFrame, align: dict | None = None) -> str:
    align = align or {}
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---:" if align.get(c) == "r" else "---" for c in cols) + "|"]
    for row in df.itertuples(index=False):
        out.append("| " + " | ".join(str(v).replace("|", "/") for v in row) + " |")
    return "\n".join(out)


def _f(x, d=3):
    return "—" if (x is None or (isinstance(x, float) and np.isnan(x))) else f"{x:.{d}f}"


def _p(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    return f"{x:.2g}" if x >= 1e-4 else f"{x:.1e}"


def other_target_tables() -> dict[str, str]:
    """Tables for the POST HOC other-target null (D54)."""
    d = pd.read_csv(OUT / "other_target_null.csv")
    sm = json.loads((OUT / "other_target_null_summary.json").read_text())
    ns = d[d["network_specific_nmet"]].copy()
    ns["verdict"] = np.where(ns["too_few_to_judge"], "too few to judge", np.where(ns["exceeds_p95"], "**above 95th**", "not above"))
    ns["_o"] = ns["too_few_to_judge"].astype(int)
    ns = ns.sort_values(["_o", "dr2_m1"], ascending=[True, False])
    tab = pd.DataFrame({
        "Target": ns["feature"], "Group": ns["group"], "n sets": ns["n_eligible"], "ΔR² (own set)": ns["dr2_m1"].map(_f),
        "Other sets: median": ns["other_median"].map(_f), "95th pct": ns["other_p95"].map(_f), "best": ns["other_max"].map(_f),
        "Rank (pct)": ns["percentile_rank"].map(lambda x: f"{x:.0f}"), "Empirical p": ns["emp_p"].map(_p), "Verdict": ns["verdict"]})
    groups = pd.DataFrame([{"Group": g, "Network-specific": v["n_network_specific"], "Judged (≥ 10 sets)": v["n_judged"],
                            "Above 95th percentile": v["n_exceed_p95"], "Too few to judge": v["n_too_few"],
                            "Median rank (pct)": "—" if v["median_percentile_rank"] is None else f"{v['median_percentile_rank']:.0f}"}
                           for g, v in sm["by_group"].items()])
    tot = {"Group": "**All**", "Network-specific": sm["n_network_specific"], "Judged (≥ 10 sets)": sm["network_specific_judged"],
           "Above 95th percentile": sm["network_specific_exceed_p95"], "Too few to judge": sm["n_targets_too_few"] - sum(
               1 for _ in [] ), "Median rank (pct)": f"{sm['median_percentile_rank_network_specific']:.0f}"}
    tot["Too few to judge"] = len(sm["network_specific_too_few"])
    groups = pd.concat([groups, pd.DataFrame([tot])], ignore_index=True)
    toofew = d[d["too_few_to_judge"]].sort_values("n_eligible")
    few = pd.DataFrame({"Target": toofew["feature"], "Network-specific (primary)": toofew["network_specific_nmet"].map(lambda v: "yes" if v else "no"),
                        "n sets": toofew["n_eligible"], "ΔR² (own set)": toofew["dr2_m1"].map(_f), "Best other set": toofew["other_max"].map(_f)})
    thin = pd.DataFrame([{"Target": r["feature"], "n sets": r["n_eligible"], "ΔR² (own set)": _f(r["dr2_m1"]),
                          "Best other set": _f(r["best_comparator"]), "Empirical p": _p(r["emp_p"])} for r in sm["resolution"]["below_best_comparator"]])
    return {"OTHER_GROUPS": _md(groups), "OTHER_TABLE": _md(tab, {"n sets": "r"}), "OTHER_TOOFEW": _md(few), "OTHER_THIN": _md(thin)}


def tables() -> dict[str, str]:
    t = pd.read_csv(OUT / "per_target_main.csv")
    g = json.loads((OUT / "global_tests.json").read_text())
    tg = t[t["role"] == "target"].sort_values("dr2_m1", ascending=False).reset_index(drop=True)

    per = pd.DataFrame({
        "#": range(1, len(tg) + 1), "Target": tg["feature"], "Class": tg["class"].str.replace("medium-supplied", "medium"),
        "k": tg["k"].astype(int), "R²(M0)": tg["r2_m0"].map(_f), "ΔR²(M1)": tg["dr2_m1"].map(_f),
        "repeats min–max": [f"{a:.3f}–{b:.3f}" for a, b in zip(tg["dr2_m1_min_rep"], tg["dr2_m1_max_rep"])],
        "ΔR²(M2)": tg["dr2_m2"].map(_f), "capture": tg["capture_ratio"].map(lambda x: _f(x, 2)),
        "N-met median / p95": [f"{a:.3f} / {b:.3f}" for a, b in zip(tg["nmet_median"], tg["nmet_p95"])],
        "p N-met": tg["p_nmet"].map(_p), "q N-met": tg["q_nmet"].map(_p),
        "N-all median": tg["nall_median"].map(_f), "p N-all": tg["p_nall"].map(_p), "q N-all": tg["q_nall"].map(_p),
        "ρ(M1)": tg["rho_m1"].map(lambda x: _f(x, 2)), "P": tg["dr2_perm"].map(lambda x: _f(x, 3)),
        "net-spec": tg["network_specific_nmet"].map(lambda v: "**yes**" if v else ""),
        "ctrl": tg["positive_control"].map(lambda v: "H5" if v else ""),
    })

    def block(name, b):
        rows = [("Targets", b["n_targets"]), ("Wilcoxon p, ΔR²(M1) − median null > 0", _p(b["wilcoxon_p_vs_median_null"])),
                ("Median excess over null", _f(b["median_excess_over_null"])),
                ("Targets above own 95th percentile (expected 5%)", f"{b['n_above_own_p95']} ({b['expected_above_p95']:.1f} expected); binomial p = {_p(b['binomial_p'])}"),
                ("Network-specific (q < 0.10 and ΔR² ≥ 0.02)", f"{b['n_network_specific']}: {', '.join(b['network_specific_targets']) or '—'}"),
                ("Classification", f"**{b['classification']}**")]
        return _md(pd.DataFrame(rows, columns=[name, "Value"]))

    h1 = g["H1"]
    H1 = _md(pd.DataFrame([("Wilcoxon p, ΔR²(M1) > 0", _p(h1["wilcoxon_p"])), ("Median ΔR²(M1)", _f(h1["median_dr2_m1"])),
                           ("Targets with ΔR²(M1) ≥ 0.02", f"{h1['n_dr2_ge_0.02']} / {g['n_targets']} ({h1['frac_dr2_ge_0.02']:.0%}; criterion ≥ 25%)"),
                           ("Targets with ΔR²(M1) > 0", h1["n_dr2_positive"]), ("H1", "**supported**" if h1["supported"] else "**not supported**")],
                          columns=["H1 — network genes add information beyond lineage + medium", "Value"]))
    H2 = block("H2 (primary) — vs size-matched random metabolic genes (N-met)", g["H2"])
    H3 = block("H3 — vs size-matched random genes from the whole universe (N-all)", g["H3"])
    c = g["H2_clustered_secondary"]
    H2C = _md(pd.DataFrame([("Clusters (single linkage, Jaccard > 0.5)", f"{c['n_clusters']} ({c['n_multi_target_clusters']} with > 1 target)"),
                            ("Wilcoxon p on cluster medians of ΔR²(M1) − median N-met", _p(c["wilcoxon_p"])),
                            ("Median cluster excess", _f(c["median_excess"])),
                            ("Clusters above own 95th N-met percentile", f"{c['n_clusters_above_p95']} / {c['n_clusters']}; binomial p = {_p(c['binomial_p'])}")],
                           columns=["H2 secondary clustered reading (D37)", "Value"]))
    h4 = g["H4"]
    H4 = _md(pd.DataFrame([("Paired Wilcoxon p, ΔR²(M2) − ΔR²(M1) > 0", _p(h4["wilcoxon_p_m2_gt_m1"])), ("Median ΔR²(M2)", _f(h4["median_dr2_m2"])),
                           ("Targets with ΔR²(M2) ≥ ΔR²(M1)", f"{h4['n_m2_ge_m1']} / {g['n_targets']}"),
                           ("Targets where M1 beats M2 by > 0.02 (scrutinise)", f"{h4['n_m1_beats_m2_by_0.02']}: {', '.join(h4['m1_beats_m2_targets']) or '—'}"),
                           ("Capture ratio ΔR²(M1)/ΔR²(M2), where ΔR²(M2) > 0.02", f"defined for {h4['n_capture_defined']} targets; median {_f(h4['median_capture_ratio'], 2)}"),
                           ("Targets with capture ratio ≥ 0.7", f"{h4['n_capture_ge_0.7']}: {', '.join(h4['capture_ge_0.7_targets']) or '—'}"),
                           ("H4", "**supported**" if h4["supported"] else "**not supported**")],
                          columns=["H4 — all-gene model at least as good as the network model", "Value"]))
    h5 = g["H5"]
    un = {r["feature"]: r for r in h5.get("unadjusted", [])}
    H5 = _md(pd.DataFrame([{"Control": r["feature"], "ΔR²(M1) adjusted": _f(r["dr2_m1"]), "p vs N-met": _p(r["p_nmet"]), "q": _p(r["q_nmet"]),
                            "R² genes only (unadjusted)": _f(un.get(r["feature"], {}).get("r2_genes_only_unadjusted", np.nan)),
                            "R²(M0)": _f(un.get(r["feature"], {}).get("r2_m0", np.nan)),
                            "Passes (ΔR² ≥ 0.02 and p < 0.05)": "**yes**" if r["passes"] else "no"} for r in h5["controls"]]))
    h6 = g["H6"]
    H6 = _md(pd.DataFrame([{"Class": c, "n": v["n"], "median R²(M0)": _f(v["median_r2_m0"]), "median ΔR²(M1)": _f(v["median_dr2_m1"]),
                            "median excess over N-met": _f(v["median_excess_over_nmet"]), "network-specific": v["n_network_specific"]}
                           for c, v in h6["by_class"].items()]))
    diag = g["pipeline_diagnostic"]
    DIAG = _md(pd.DataFrame([{"Diagnostic": d["feature"], "k": d["k"], "ΔR²(M1)": _f(d["dr2_m1"]), "N-met median": _f(d["nmet_median"]),
                              "p vs N-met": _p(d["p_nmet"]), "Status": d["note"]} for d in diag]))
    s = g["sensitivity"]
    SENS = _md(pd.DataFrame([{"Analysis": v["analysis"], "n targets": v["n_targets"], "Wilcoxon p": _p(v["wilcoxon_p_vs_median_null"]),
                              "median excess": _f(v["median_excess_over_null"]), "above p95": v["n_above_own_p95"],
                              "network-specific": f"{v['n_network_specific']}" + (f" ({', '.join(v['network_specific_targets'])})" if v.get("network_specific_targets") else ""),
                              "H2 classification": f"**{v['classification']}**"} for k, v in s.items() if k != "robustness"]))
    rb = s["robustness"]
    ROBUST = f"**{'Robust' if rb['robust'] else 'Not robust'}**: the H2 classification is *{rb['H2_main']}* in the main analysis and " + \
             ("unchanged in all " if rb["robust"] else "changes in at least one of ") + f"{rb['n_analyses_available']} sensitivity analyses."
    s1 = g["S1_permuted_control"]
    S1 = f"median ΔR²(P) = {s1['median']:+.4f}, max = {s1['max']:+.4f}, targets above 0.02: {s1['n_above_0.02']} → **{'PASS' if s1['pass'] else 'FAIL'}**"
    extra = other_target_tables() if (OUT / "other_target_null.csv").exists() else {}
    return {**extra, "PER_TARGET": _md(per, {"#": "r", "k": "r"}), "H1": H1, "H2": H2, "H2C": H2C, "H3": H3, "H4": H4, "H5": H5, "H6": H6,
            "DIAG": DIAG, "SENS": SENS, "ROBUST": ROBUST, "S1": S1}


def main() -> None:
    tpl = (REPO / "docs" / "phase3_report.template.md").read_text()
    for k, v in tables().items():
        tpl = tpl.replace("{{" + k + "}}", v)
    assert "{{" not in tpl, "unfilled placeholder"
    (REPO / "docs" / "phase3_report.md").write_text(tpl)
    print("wrote docs/phase3_report.md")


if __name__ == "__main__":
    main()
