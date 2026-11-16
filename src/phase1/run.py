"""Phase 1: network, gene mapping, metabolite mapping, gene sets and outcome-blind target selection.

    .venv/bin/python -m src.phase1.run

Nothing here uses metabolite-expression relationships: the only data inputs to target selection are
the mapping (IDs), the Human-GEM network, expression *levels* (the expressed-gene filter) and each
metabolite's own pre-imputation missingness and spread.
"""
from __future__ import annotations

import json
from itertools import combinations

import numpy as np
import pandas as pd

from ..phase0.download import hash_file
from ..phase0.sources import HUMANGEM_COMMIT, HUMANGEM_TAG, PROCESSED, REPO
from . import genemap, genesets, humangem, metmap

OUT = PROCESSED / "phase1"
CURATION = REPO / "data" / "curation"
K_MIN, K_MAX = 3, 60                 # research_plan §5 rule 3 (D4)
MAX_MISSING = 0.20                   # rule 4 (D18)
MIN_SD = 0.05                        # rule 4 "non-trivial variance" (D25), log10 units

POSITIVE_CONTROLS = {                # expectations H5 as amended 2026-11-16 (D34, option b)
    # feature: (literature genes, in-set genes named in the amended H5)
    "1-methylnicotinamide": (["NNMT"], ["NNMT"]),
    "kynurenine": (["IDO1", "TDO2"], ["AFMID", "KYNU"]),
    "creatine": (["SLC6A8", "GATM"], ["GAMT", "SLC6A8"]),
    "N-carbamoyl-beta-alanine": (["DPYD"], ["DPYS", "UPB1"]),
    "taurine": (["SLC6A6"], ["SLC6A6"]),
}
PIPELINE_DIAGNOSTICS = ["1-methylnicotinamide"]   # D34: M1 fitted outside the selection rules; not a target
JACCARD_CLUSTER = 0.5                              # D37: single-linkage clusters for the secondary global reading
K_MIN_ENZYMATIC = 3                                # D36: sensitivity analysis 2 drops targets below this


def _write(df: pd.DataFrame, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / name, index=False, float_format="%.6g", lineterminator="\n")


def run(write: bool = True) -> dict:
    cells = pd.read_csv(PROCESSED / "cell_lines_phase0.csv")
    primary = cells.loc[cells["in_primary_set"], "CCLE_ID"].tolist()
    inv = pd.read_csv(PROCESSED / "metabolite_inventory.csv")
    polar = inv[inv["assay_class"].str.startswith("polar")].set_index("metabolite")

    # ---- metabolite mapping
    mapping, species = metmap.build()
    mapping["single_hg_target"] = mapping["hg_target_proposed"]

    # ---- gene mapping and expressed genes
    gtf = genemap.read_gtf_genes()
    eg, _ = genemap.expression_genes(primary)
    gmap = genemap.map_genes(eg, gtf)
    expressed_ids = set(eg.loc[eg["expressed"], "gene_id"])
    hg2ccle = dict(zip(gmap.loc[gmap["mapped"], "hg_gene"], gmap.loc[gmap["mapped"], "ccle_gene_id"]))
    sym = dict(zip(gmap["hg_gene"], gmap["hg_symbol"]))

    # ---- gene sets for every feature with a single Human-GEM target (confident or provisional)
    model = humangem.load_model()
    hg_names = {m.id[:-1]: m.name for m in model.metabolites}
    gs_rows, members = [], []
    sets: dict[str, dict[str, set]] = {}
    for r in mapping[mapping["single_hg_target"] != ""].itertuples():
        g = genesets.gene_set(model, r.single_hg_target)
        rec = {"feature": r.feature, "hg_metsNoComp": g.target, "hg_name": hg_names.get(g.target, ""),
               "compartments": g.compartments, "n_reactions_all": g.n_reactions_all,
               "n_excluded": g.n_excluded, "excluded_by": json.dumps(g.excluded_by, sort_keys=True),
               "n_retained": g.n_retained, "n_transport": g.n_transport, "n_enzymatic": g.n_enzymatic,
               "n_retained_no_gpr": g.n_retained_no_gpr, "n_enzymatic_no_gpr": g.n_enzymatic_no_gpr}
        sets[r.feature] = {}
        for ver, genes in (("primary", g.genes_primary), ("enzymatic", g.genes_enzymatic)):
            mapped = {hg2ccle[x] for x in genes if x in hg2ccle}
            expd = mapped & expressed_ids
            rec[f"k_hg_{ver}"] = len(genes)
            rec[f"k_mapped_{ver}"] = len(mapped)
            rec[f"k_expressed_{ver}"] = len(expd)
            rec[f"frac_lost_unmapped_{ver}"] = (1 - len(mapped) / len(genes)) if genes else np.nan
            sets[r.feature][ver] = expd
            for x in sorted(genes):
                members.append({"feature": r.feature, "version": ver, "hg_gene": x, "symbol": sym.get(x, ""),
                                "ccle_gene_id": hg2ccle.get(x, ""), "mapped": x in hg2ccle,
                                "expressed": hg2ccle.get(x, "") in expressed_ids})
        gs_rows.append(rec)
    gs = pd.DataFrame(gs_rows)
    members = pd.DataFrame(members)

    # ---- target rules with attrition (sequential)
    t = mapping[["feature", "status", "needs_review", "single_hg_target"]].copy()
    t = t.merge(gs, on="feature", how="left")
    t["raw_missing_frac_primary"] = t["feature"].map(polar["raw_missing_frac_primary"])
    t["sd_primary"] = t["feature"].map(polar["sd_primary"])
    t["rule1_hilic"] = t["feature"].isin(polar.index)
    t["rule2_confident_strict"] = t["status"].eq("1:1 confident")
    t["rule2_confident_incl_provisional"] = t["status"].isin(["1:1 confident", "1:1 confident (provisional)"]) \
        & (t["single_hg_target"] != "")
    t["rule3_size"] = t["k_expressed_primary"].between(K_MIN, K_MAX)
    t["rule4_missing"] = t["raw_missing_frac_primary"] <= MAX_MISSING
    t["rule4_variance"] = t["sd_primary"] >= MIN_SD
    t["target_strict"] = t[["rule1_hilic", "rule2_confident_strict", "rule3_size", "rule4_missing",
                            "rule4_variance"]].all(axis=1)
    t["target_if_review_accepted"] = t[["rule1_hilic", "rule2_confident_incl_provisional", "rule3_size",
                                        "rule4_missing", "rule4_variance"]].all(axis=1)
    t["depends_on_review"] = t["target_if_review_accepted"] & ~t["target_strict"]

    def attrition(col2):
        steps, m = [], pd.Series(True, index=t.index)
        steps.append(("CCLE metabolomics features (all assays)", len(inv)))
        for label, col in (("Rule 1: measured by HILIC (polar)", "rule1_hilic"),
                           ("Rule 2: 1:1 confident mapping", col2),
                           (f"Rule 3: {K_MIN}-{K_MAX} expressed genes (primary gene set)", "rule3_size"),
                           (f"Rule 4a: <= {MAX_MISSING:.0%} missing before imputation", "rule4_missing"),
                           (f"Rule 4b: SD >= {MIN_SD} log10 units", "rule4_variance")):
            m &= t[col].fillna(False)
            steps.append((label, int(m.sum())))
        return steps
    att = pd.DataFrame({"step": [s for s, _ in attrition("rule2_confident_strict")],
                        "n_strict": [n for _, n in attrition("rule2_confident_strict")],
                        "n_if_review_accepted": [n for _, n in attrition("rule2_confident_incl_provisional")]})

    # marginal failure reasons for mapped single-target features
    def fail_reason(r):
        if not r.rule2_confident_incl_provisional:
            return f"mapping: {r.status}"
        out = []
        if not r.rule3_size:
            k = r.k_expressed_primary
            out.append(f"gene set {int(k)} expressed genes ({'< 3' if k < K_MIN else '> 60'})")
        if not r.rule4_missing:
            out.append(f"{r.raw_missing_frac_primary:.0%} missing")
        if not r.rule4_variance:
            out.append(f"SD {r.sd_primary:.3f}")
        return "; ".join(out) if out else ("passes (needs mapping review)" if r.depends_on_review else "passes")
    t["selection_result"] = t.apply(fail_reason, axis=1)

    # ---- annotations for candidate targets (curated, outcome-blind)
    ann = pd.read_csv(CURATION / "target_annotations.tsv", sep="\t", dtype=str).fillna("")
    t = t.merge(ann, on="feature", how="left")
    t["positive_control"] = t["feature"].isin(POSITIVE_CONTROLS)
    t["pipeline_diagnostic"] = t["feature"].isin(PIPELINE_DIAGNOSTICS) & ~t["target_strict"]
    t["in_sensitivity2_enzymatic"] = t["target_strict"] & (t["k_expressed_enzymatic"] >= K_MIN_ENZYMATIC)

    def _genes_status(f, genes):
        mp = members[(members["feature"] == f) & (members["version"] == "primary")]
        return ", ".join(f"{g}:{'expressed' if ((mp.symbol == g) & mp.expressed).any() else ('in set, not expressed' if (mp.symbol == g).any() else 'NOT in set')}"
                         for g in genes)
    pc_rows = []
    for f, (lit, inset) in POSITIVE_CONTROLS.items():
        r = t[t["feature"] == f].iloc[0]
        pc_rows.append({"feature": f, "literature_genes": ", ".join(lit), "literature_genes_status": _genes_status(f, lit),
                        "h5_in_set_genes": ", ".join(inset), "h5_in_set_genes_status": _genes_status(f, inset),
                        "mapping_status": r["status"], "k_expressed_primary": r["k_expressed_primary"],
                        "selected": bool(r["target_strict"]), "selection_result": r["selection_result"],
                        "role": ("target + H5 control" if r["target_strict"] else
                                 ("pipeline diagnostic (not a target)" if r["pipeline_diagnostic"] else "excluded by rule 3"))})
    pcs = pd.DataFrame(pc_rows)
    diag = members[members["feature"].isin(t.loc[t["pipeline_diagnostic"], "feature"]) & (members["version"] == "primary")]

    # ---- overlap between targets (primary, expressed genes)
    tg = t.loc[t["target_strict"], "feature"].tolist()
    jac = pd.DataFrame(np.eye(len(tg)), index=tg, columns=tg)
    pairs = []
    for a, b in combinations(tg, 2):
        j = genesets.jaccard(sets[a]["primary"], sets[b]["primary"])
        jac.loc[a, b] = jac.loc[b, a] = j
        pairs.append({"a": a, "b": b, "jaccard": j,
                      "n_shared": len(sets[a]["primary"] & sets[b]["primary"])})
    pairs = pd.DataFrame(pairs).sort_values("jaccard", ascending=False)

    # D37: single-linkage clusters at Jaccard > 0.5 (= connected components of the > 0.5 graph)
    parent = {f: f for f in tg}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b, j in zip(pairs["a"], pairs["b"], pairs["jaccard"]):
        if j > JACCARD_CLUSTER:
            parent[find(a)] = find(b)
    roots = {f: find(f) for f in tg}
    comp: dict[str, list] = {}
    for f, rt in roots.items():
        comp.setdefault(rt, []).append(f)
    comps = sorted(comp.values(), key=lambda c: (-len(c), sorted(c)[0]))
    cid = {f: f"C{i + 1:02d}" for i, c in enumerate(comps) for f in c}
    clusters = pd.DataFrame({"feature": tg, "cluster_id": [cid[f] for f in tg],
                             "cluster_size": [len(comp[roots[f]]) for f in tg]}).sort_values(["cluster_id", "feature"])

    # D33 report: genes that a single subsystem contributes to > 10 targets
    contrib: dict[tuple[str, str], set] = {}
    for f in tg:
        g = genesets.gene_set(model, mapping.loc[mapping["feature"] == f, "hg_target_proposed"].iloc[0])
        for rid in g.reactions_retained:
            rx = model.reactions.get_by_id(rid)
            for gene in rx.genes:
                if hg2ccle.get(gene.id) in expressed_ids:
                    contrib.setdefault((rx.subsystem, gene.id), set()).add(f)
    shared = pd.DataFrame([{"subsystem": sub, "hg_gene": gid, "symbol": sym.get(gid, ""), "n_targets": len(fs),
                            "targets": "; ".join(sorted(fs))} for (sub, gid), fs in contrib.items() if len(fs) > 10])
    shared = shared.sort_values(["subsystem", "n_targets", "symbol"], ascending=[True, False, True]) if len(shared) else shared

    summary = {
        "humangem": {"tag": HUMANGEM_TAG, "commit": HUMANGEM_COMMIT, "n_reactions": len(model.reactions),
                     "n_metabolites": len(model.metabolites), "n_genes": len(model.genes)},
        "gene_mapping": {"n_hg_genes": int(len(gmap)), "n_mapped": int(gmap["mapped"].sum()),
                         "by_ensembl": int((gmap["method"] == "ensembl").sum()),
                         "by_symbol": int((gmap["method"] == "symbol").sum()),
                         "n_unmapped": int((~gmap["mapped"]).sum()),
                         "mapping_rate": float(gmap["mapped"].mean()),
                         "n_mapped_expressed_primary": int(gmap["ccle_gene_id"].isin(expressed_ids).sum()),
                         "n_expressed_genes_primary": int(len(expressed_ids))},
        "metabolite_mapping": {k: int(v) for k, v in mapping["status"].value_counts().items()} |
                              {"n_needs_review": int(mapping["needs_review"].sum())},
        "targets": {"n_targets": int(t["target_strict"].sum()),
                    "n_if_review_accepted": int(t["target_if_review_accepted"].sum()),
                    "gate_min_targets": 15, "gate_pass": bool(t["target_strict"].sum() >= 15),
                    "n_pipeline_diagnostics": int(t["pipeline_diagnostic"].sum()),
                    "n_in_sensitivity2_enzymatic": int(t["in_sensitivity2_enzymatic"].sum()),
                    "dropped_from_sensitivity2": sorted(t.loc[t["target_strict"] & ~t["in_sensitivity2_enzymatic"], "feature"]),
                    "n_clusters_jaccard_gt_0.5": int(clusters["cluster_id"].nunique()),
                    "n_multi_target_clusters": int((clusters.groupby("cluster_id").size() > 1).sum()),
                    "class_counts": {k: int(v) for k, v in t.loc[t["target_strict"], "class"].value_counts().items()}},
        "overlap": {"n_pairs": int(len(pairs)),
                    "median_jaccard": float(pairs["jaccard"].median()) if len(pairs) else None,
                    "n_pairs_jaccard_gt_0.5": int((pairs["jaccard"] > 0.5).sum()) if len(pairs) else 0,
                    "n_pairs_any_shared": int((pairs["n_shared"] > 0).sum()) if len(pairs) else 0},
    }
    out = {"mapping": mapping, "species": species, "gene_mapping": gmap, "genesets": gs, "members": members,
           "targets": t, "attrition": att, "positive_controls": pcs, "jaccard": jac, "jaccard_pairs": pairs,
           "clusters": clusters, "shared_subsystem_genes": shared, "diagnostic_members": diag, "summary": summary}
    if write:
        from . import idlookup
        idlookup.register_manifest()
        _write(mapping, "metabolite_mapping.csv")
        rec = pd.read_csv(CURATION / "review_recommendations.tsv", sep="\t", dtype=str)
        rv = mapping[mapping["review_decision"] != ""][
            ["feature", "status_rule_based", "reason", "review_decision", "status", "review_rationale", "review_ref"]]
        rv = rv.merge(rec[["feature", "recommendation"]], on="feature", how="left").merge(
            t[["feature", "target_strict", "selection_result"]], on="feature", how="left")
        rv = rv.rename(columns={"status": "status_final", "reason": "issue"})
        assert not mapping["needs_review"].any(), "unresolved review items remain"
        _write(rv, "metabolite_mapping_review_decisions.csv")
        (OUT / "metabolite_mapping_needs_review.csv").unlink(missing_ok=True)
        _write(gmap, "gene_mapping.csv")
        _write(gmap[~gmap["mapped"]], "gene_mapping_unmapped.csv")
        _write(gs, "gene_set_summary.csv")
        _write(members, "gene_set_members.csv")
        _write(t, "target_selection.csv")
        _write(att, "target_attrition.csv")
        _write(pcs, "positive_controls.csv")
        _write(diag, "pipeline_diagnostic_gene_sets.csv")
        _write(clusters, "target_clusters.csv")
        _write(shared, "shared_subsystem_genes.csv")
        jac.round(4).to_csv(OUT / "target_jaccard_matrix.csv", lineterminator="\n")
        _write(pairs, "target_jaccard_pairs.csv")
        (OUT / "phase1_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        rows = [{"file": str(p.relative_to(REPO)), "sha256": hash_file(p)[0]} for p in sorted(OUT.glob("*"))
                if p.name != "FREEZE.tsv"]
        pd.DataFrame(rows).to_csv(OUT / "FREEZE.tsv", sep="\t", index=False, lineterminator="\n")
    return out


if __name__ == "__main__":
    o = run()
    print(json.dumps(o["summary"], indent=2))
