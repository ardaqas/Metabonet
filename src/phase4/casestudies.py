"""Phase 4 case studies (D55): gene roles, ridge coefficients and out-of-fold predictions for three targets.

    .venv/bin/python -m src.phase4.casestudies

Cases (fixed by rule before any coefficient was computed): top target (2-deoxycytidine), best passing positive
control (creatine), the failing positive control (kynurenine). Everything here is descriptive and associative.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from ..phase0.sources import PROCESSED, REPO
from ..phase1 import genesets, humangem
from ..phase2 import folds as fo
from ..phase2 import pipeline as pl
from ..phase2 import ridge
from ..phase2.data import baseline_design, load_inputs
from ..phase2.run import FOLDS

OUT = REPO / "results" / "phase4"
CASES = {"2-deoxycytidine": "top target", "creatine": "positive control (passes)", "kynurenine": "positive control (fails)"}


def gene_roles(model, nocomp: str, retained: list[str]) -> dict[str, dict]:
    """Per Human-GEM gene: role(s) relative to the target over the retained reactions of its first-neighbour set.

    producing  = target is a product in the reaction's forward direction and the reaction is irreversible
    consuming  = target is a reactant in the forward direction and the reaction is irreversible
    reversible = reaction can run both ways (direction not identifiable from expression; research_plan §3.3)
    transport  = target on both sides after collapsing compartments (D24)
    """
    roles: dict[str, dict] = {}
    for rid in retained:
        rx = model.reactions.get_by_id(rid)
        if not rx.genes:
            continue
        lhs = {genesets.base_id(m.id) for m, c in rx.metabolites.items() if c < 0}
        rhs = {genesets.base_id(m.id) for m, c in rx.metabolites.items() if c > 0}
        if nocomp in lhs and nocomp in rhs:
            role = "transport"
        elif rx.lower_bound < 0 < rx.upper_bound:
            role = "reversible"
        elif rx.upper_bound > 0:
            role = "producing" if nocomp in rhs else "consuming"
        else:   # runs only backwards
            role = "consuming" if nocomp in rhs else "producing"
        label = (rx.name or rid).strip() or rid
        for g in rx.genes:
            d = roles.setdefault(g.id, {"roles": set(), "reactions": []})
            d["roles"].add(role)
            d["reactions"].append(f"{label} [{rx.subsystem}]")
    return roles


def fit_all_folds(ev: pl.Evaluator, y: np.ndarray, tf: pl.TargetFolds, gidx: np.ndarray):
    """M1 refitted on every training fold exactly as in the pipeline; returns coefficients and repeat-0 OOF predictions."""
    coefs, alphas = [], []
    pred1 = np.empty(ev.n)
    pred0 = np.empty(ev.n)
    r2 = []
    sst = ((y - y.mean()) ** 2).sum()
    for r, reps in enumerate(ev.folds):
        pr = np.empty(ev.n)
        for f, (tr, te) in enumerate(reps):
            Ztr, Zte = ev._z(r, f, tr, gidx), ev._z(r, f, te, gidx)
            coef, b0, a, _ = ridge.fit(Ztr, tf.res_tr[r][f], ev.alphas)
            coefs.append(coef)
            alphas.append(a)
            pr[te] = tf.base_te[r][f] + ridge.predict(Zte, coef, b0)
            if r == 0:
                pred0[te] = tf.base_te[r][f]
        if r == 0:
            pred1 = pr.copy()
        r2.append(1 - ((y - pr) ** 2).sum() / sst)
    return np.array(coefs), np.array(alphas), pred0, pred1, np.array(r2)


def full_data_fit(inp, y: np.ndarray, D: np.ndarray, gidx: np.ndarray):
    """One descriptive fit on all 735 lines: OLS baseline, genes standardised on all lines, ridge by LOO."""
    beta, *_ = np.linalg.lstsq(D, y, rcond=None)
    res = y - D @ beta
    X = inp.X[:, gidx].astype(np.float64)
    sd = X.std(axis=0); sd[sd == 0] = 1.0
    Z = (X - X.mean(axis=0)) / sd
    coef, _, a, _ = ridge.fit(Z, res)
    return coef, float(ridge.ALPHAS[a])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inp = load_inputs("primary")
    ev = pl.Evaluator(inp, fo.fold_indices(pd.read_csv(FOLDS), inp.cells))
    D, _ = baseline_design(inp.cells)
    model = humangem.load_model()
    sel = pd.read_csv(PROCESSED / "phase1" / "target_selection.csv").set_index("feature")
    members = pd.read_csv(PROCESSED / "phase1" / "gene_set_members.csv")
    gmap = pd.read_csv(PROCESSED / "phase1" / "gene_mapping.csv").set_index("ccle_gene_id")
    main_df = pd.read_csv(REPO / "results" / "phase3" / "per_target_main.csv").set_index("feature")
    summary = {}
    for name, why in CASES.items():
        y = inp.Y[name].to_numpy(dtype=float)
        tf = ev.target_folds(y)
        gene_ids = inp.sets_primary[name]
        gidx = inp.idx(gene_ids)
        coefs, alphas, pred0, pred1, r2 = fit_all_folds(ev, y, tf, gidx)
        z = np.load(REPO / "results" / "phase3" / "main" / f"{name}.npz")
        assert np.allclose(r2, z["r2_m1"], atol=1e-12), f"{name}: refit does not reproduce Phase 3 R2"
        coef_full, alpha_full = full_data_fit(inp, y, D, gidx)

        nocomp = sel.at[name, "single_hg_target"]
        retained = genesets.gene_set(model, nocomp).reactions_retained
        roles = gene_roles(model, nocomp, retained)
        mem = members[(members["feature"] == name) & (members["version"] == "primary") & members["expressed"]]
        rows = []
        for gi, gid in enumerate(gene_ids):
            hg = mem.loc[mem["ccle_gene_id"] == gid, "hg_gene"].iloc[0]
            sym = mem.loc[mem["ccle_gene_id"] == gid, "symbol"].iloc[0]
            col = coefs[:, gi]
            expr = inp.X[:, gidx[gi]].astype(float)
            r = roles.get(hg, {"roles": set(), "reactions": []})
            rows.append({"feature": name, "gene": sym, "hg_gene": hg, "ccle_gene_id": gid,
                         "roles": "/".join(sorted(r["roles"])) or "no GPR role found", "n_reactions": len(r["reactions"]),
                         "example_reactions": " ; ".join(sorted(set(r["reactions"]))[:3]),
                         "coef_mean": col.mean(), "coef_sd": col.std(ddof=1), "sign_consistency": max((col > 0).mean(), (col < 0).mean()),
                         "coef_full_data": coef_full[gi], "expr_mean_log2tpm1": expr.mean(), "expr_sd_log2tpm1": expr.std(ddof=1)})
        gdf = pd.DataFrame(rows)
        gdf["abs_coef_mean"] = gdf["coef_mean"].abs()
        gdf = gdf.sort_values("abs_coef_mean", ascending=False).drop(columns="abs_coef_mean").reset_index(drop=True)
        gdf.to_csv(OUT / f"case_{name}_genes.csv", index=False, float_format="%.5g", lineterminator="\n")
        pd.DataFrame({"CCLE_ID": inp.cells["CCLE_ID"], "lineage": inp.cells["lineage_collapsed"], "observed": y,
                      "pred_m0": pred0, "pred_m1": pred1}).to_csv(OUT / f"case_{name}_oof_repeat0.csv", index=False,
                                                                  float_format="%.5g", lineterminator="\n")
        summary[name] = {"why": why, "k": len(gene_ids), "hg_metabolite": nocomp, "n_retained_reactions": len(retained),
                         "r2_m0": float(tf.r2_m0.mean()), "dr2_m1": float((r2 - tf.r2_m0).mean()), "dr2_m2": float(main_df.at[name, "dr2_m2"]),
                         "p_nmet": float(main_df.at[name, "p_nmet"]), "alpha_median_idx": int(np.median(alphas)),
                         "alpha_full_data": alpha_full, "reproduces_phase3_r2": True,
                         "oof_repeat0_pearson_m1": float(np.corrcoef(y, pred1)[0, 1]), "oof_repeat0_pearson_m0": float(np.corrcoef(y, pred0)[0, 1])}
        print(f"{name}: k={len(gene_ids)} dR2={summary[name]['dr2_m1']:.4f} (reproduces Phase 3)")
    (OUT / "case_studies_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
