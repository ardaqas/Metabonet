"""Render Phase 1 markdown tables from data/processed/phase1/ into docs/phase1_report.md.

    .venv/bin/python -m src.phase1.report_tables

docs/phase1_report.template.md holds the hand-written text with {{PLACEHOLDER}} lines; this script
replaces each placeholder with a table generated from the frozen CSV outputs.
"""
from __future__ import annotations

import json

import pandas as pd

from ..phase0.sources import PROCESSED, REPO

D = PROCESSED / "phase1"


def _first(s: str) -> str:
    return str(s).split("|")[0].split(",")[0].strip()


def _md(df: pd.DataFrame, align: dict | None = None) -> str:
    align = align or {}
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join("---:" if align.get(c) == "r" else "---" for c in cols) + "|"]
    for row in df.itertuples(index=False):
        out.append("| " + " | ".join(str(v).replace("|", "/") for v in row) + " |")
    return "\n".join(out)


def tables() -> dict[str, str]:
    mp = pd.read_csv(D / "metabolite_mapping.csv", dtype=str).fillna("")
    t = pd.read_csv(D / "target_selection.csv").fillna("")
    rv = pd.read_csv(D / "metabolite_mapping_review_decisions.csv").fillna("")
    att = pd.read_csv(D / "target_attrition.csv")
    pcs = pd.read_csv(D / "positive_controls.csv").fillna("")
    gun = pd.read_csv(D / "gene_mapping_unmapped.csv").fillna("")
    cl = pd.read_csv(D / "target_clusters.csv")
    sh = pd.read_csv(D / "shared_subsystem_genes.csv").fillna("")
    s = json.loads((D / "phase1_summary.json").read_text())

    hgm = pd.read_csv(REPO / "data" / "raw" / "humangem" / "metabolites.tsv", sep="\t", dtype=str).fillna("")
    hg_hmdb = hgm.groupby("metsNoComp")["metHMDBID"].first().to_dict()
    cid = dict(zip(cl["feature"], cl["cluster_id"]))
    csize = dict(zip(cl["feature"], cl["cluster_size"]))
    tg = t[t["target_strict"].astype(bool)].merge(mp[["feature", "kegg", "hmdb", "hmdb_source"]], on="feature")
    tg["hmdb"] = [((_first(h) + (" (curator, unverified)" if "curator" in src else "")) if _first(h) else
                   (f"{hg_hmdb.get(m, '').split(';')[0]} (HG)" if hg_hmdb.get(m) else ""))
                  for h, m, src in zip(tg["hmdb"], tg["single_hg_target"], tg["hmdb_source"])]
    tg = tg.sort_values(["class", "feature"])
    targets = pd.DataFrame({
        "#": range(1, len(tg) + 1),
        "Target (CCLE feature)": tg["feature"],
        "HMDB": tg["hmdb"].replace("", "—"),
        "KEGG": tg["kegg"].map(_first).replace("", "—"),
        "Human-GEM ID": tg["single_hg_target"],
        "Human-GEM name": tg["hg_name"],
        "k (primary)": tg["k_expressed_primary"].astype(int),
        "k (enzymatic)": tg["k_expressed_enzymatic"].astype(int),
        "S2": tg["in_sensitivity2_enzymatic"].map(lambda v: "yes" if str(v) == "True" else "dropped"),
        "Class": tg["class"],
        "Pos. control": tg["positive_control"].map(lambda v: "yes" if str(v) == "True" else ""),
        "Cluster": tg["feature"].map(lambda f: cid[f] if csize.get(f, 1) > 1 else ""),
    })

    mpt = mp.merge(t[["feature", "selection_result"]], on="feature", how="left")
    mapping = pd.DataFrame({
        "Feature": mpt["feature"],
        "Status": mpt["status"] + mpt["review_decision"].map(lambda v: f" (reviewed: {v})" if v else ""),
        "KEGG": mpt["kegg"].str.replace(" | ", "; ", regex=False).replace("", "—"),
        "Human-GEM": mpt["hg_target_proposed"].where(mpt["hg_target_proposed"] != "",
                                                       mpt["hg_metsNoComp"].replace("", "—")),
        "Reason": mpt["review_rationale"].where(mpt["review_rationale"] != "", mpt["reason"]).str.slice(0, 180),
        "Selection": mpt["selection_result"].fillna("—"),
    })

    rv["_o"] = rv["review_decision"].map({"accept": 0, "exclude": 1, "keep_ambiguous": 2, "moot": 3})
    rv = rv.sort_values(["_o", "feature"])
    review = pd.DataFrame({
        "Feature": rv["feature"],
        "Rule-based status": rv["status_rule_based"],
        "Issue": rv["issue"].str.slice(0, 150),
        "Decision (D31)": rv["review_decision"],
        "Final status": rv["status_final"],
        "Rationale": rv["review_rationale"],
        "Selection": rv["selection_result"],
    })

    attrition = att.rename(columns={"step": "Step", "n_strict": "n"}).drop(columns=["n_if_review_accepted"])

    pc = pcs.rename(columns={"feature": "Control", "literature_genes_status": "Literature gene(s)",
                             "h5_in_set_genes_status": "H5 in-set gene(s)", "k_expressed_primary": "k expressed",
                             "role": "Role"})[["Control", "Literature gene(s)", "H5 in-set gene(s)", "k expressed", "Role"]]
    pc["k expressed"] = pc["k expressed"].map(lambda v: int(float(v)) if v != "" else "")

    unm = gun.rename(columns={"hg_gene": "Human-GEM gene", "hg_symbol": "Symbol", "reason": "Reason"})[
        ["Human-GEM gene", "Symbol", "Reason"]]

    status_counts = pd.DataFrame(
        [(k, v) for k, v in s["metabolite_mapping"].items() if k != "n_needs_review"],
        columns=["Status", "Features"])

    multi = cl[cl["cluster_size"] > 1].groupby("cluster_id").agg(
        n=("feature", "size"), targets=("feature", lambda x: ", ".join(sorted(x)))).reset_index()
    clusters = multi.rename(columns={"cluster_id": "Cluster", "n": "Targets (n)", "targets": "Members"})

    shared = (sh.groupby("subsystem").agg(n_genes=("symbol", "size"), max_targets=("n_targets", "max"),
                                           genes=("symbol", lambda x: ", ".join(f"{g} ({n})" for g, n in zip(x, sh.loc[x.index, "n_targets"]))))
              .reset_index().rename(columns={"subsystem": "Subsystem", "n_genes": "Genes in > 10 targets",
                                             "max_targets": "Max targets", "genes": "Genes (n targets)"})) if len(sh) else \
        pd.DataFrame({"Subsystem": ["(none)"]})

    sens2 = t[t["target_strict"].astype(bool) & ~t["in_sensitivity2_enzymatic"].astype(bool)][
        ["feature", "k_expressed_primary", "k_expressed_enzymatic"]].rename(columns={
            "feature": "Target", "k_expressed_primary": "k (primary)", "k_expressed_enzymatic": "k (enzymatic)"})
    sens2["k (primary)"] = sens2["k (primary)"].astype(int)
    sens2["k (enzymatic)"] = sens2["k (enzymatic)"].astype(int)

    return {"TARGETS": _md(targets, {"#": "r", "k (primary)": "r", "k (enzymatic)": "r"}),
            "MAPPING": _md(mapping), "REVIEW": _md(review), "ATTRITION": _md(attrition),
            "POSCONTROLS": _md(pc), "UNMAPPED_GENES": _md(unm), "STATUS_COUNTS": _md(status_counts),
            "CLUSTERS": _md(clusters), "SHARED_SUBSYSTEMS": _md(shared), "SENS2_DROPPED": _md(sens2)}


def main() -> None:
    tpl = (REPO / "docs" / "phase1_report.template.md").read_text()
    for k, v in tables().items():
        tpl = tpl.replace("{{" + k + "}}", v)
    assert "{{" not in tpl, "unfilled placeholder"
    (REPO / "docs" / "phase1_report.md").write_text(tpl)
    # expectations.md §9: frozen target table between markers (regenerated, never hand-edited)
    exp = REPO / "docs" / "expectations.md"
    txt = exp.read_text()
    a, b = "<!-- TARGETS:START -->", "<!-- TARGETS:END -->"
    assert a in txt and b in txt, "expectations.md markers missing"
    pre, rest = txt.split(a, 1)
    _, post = rest.split(b, 1)
    exp.write_text(pre + a + "\n" + tables()["TARGETS"] + "\n" + b + post)
    print("wrote docs/phase1_report.md and expectations.md §9")


if __name__ == "__main__":
    main()
