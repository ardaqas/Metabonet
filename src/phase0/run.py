"""Run the full Phase 0 audit and write frozen outputs to data/processed/.

    .venv/bin/python -m src.phase0.run

Outputs are deterministic (sorted rows, fixed float format) so that
data/processed/FREEZE.tsv hashes are reproducible.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import audit, load, match, medium
from .download import hash_file
from .sources import PROCESSED, REPO

MIN_PRIMARY_LINES = 400
MIN_MEDIUM_COVERAGE = 0.90

CELL_COLUMNS = [
    "CCLE_ID", "DepMap_ID", "DepMap_ID_metab", "DepMap_ID_anno2019", "ModelID", "CCLEName", "CellLineName", "RRID",
    "cellosaurus_ac", "cellosaurus_name", "PatientID",
    "lineage", "lineage_collapsed", "OncotreePrimaryDisease", "OncotreeSubtype", "site_primary_2019",
    "li2019_classification", "growth_pattern",
    "is_haematopoietic", "is_non_cancerous", "is_solid_tumour", "is_duplicate_line",
    "is_problematic", "problematic_type", "problematic_text",
    "donor_group_id", "donor_group_size", "donor_group_size_primary", "donor_link_reason",
    "medium_raw_li2019", "medium_raw_anno2019", "medium_base_li2019", "medium_base_anno2019",
    "medium_sources_agree", "medium_components", "medium_base", "medium_category", "medium_source",
    "serum_pct", "medium_assigned",
    "has_metabolomics", "has_rnaseq", "in_model", "in_primary_set",
]


def _write_csv(df: pd.DataFrame, name: str) -> None:
    df.to_csv(PROCESSED / name, index=False, float_format="%.6g", lineterminator="\n")


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return round(float(o), 8)
    if isinstance(o, float):
        return round(o, 8)
    if isinstance(o, tuple):
        return [_jsonable(x) for x in o]
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_jsonable(x) for x in o]
    return o


def run(write: bool = True) -> dict:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    res = match.build()
    cells = res["cells"].sort_values("CCLE_ID").reset_index(drop=True)
    issues = res["issues"].sort_values(["issue", "CCLE_ID"], na_position="first").reset_index(drop=True)
    attrition = res["attrition"]

    # medium normalisation table (all metabolomics lines) + final collapsed category per base
    norm = medium.normalisation_table(cells["medium_raw_li2019"], cells["medium_raw_anno2019"])
    base_to_cat = (cells.loc[cells["in_primary_set"]].groupby("medium_base")["medium_category"].first().to_dict())
    norm["medium_category_primary"] = norm["medium_base"].map(base_to_cat).fillna(
        norm["medium_base"].map(lambda b: "unknown" if b in medium.UNASSIGNED else "(not in primary set)"))
    prim = cells[cells["in_primary_set"]]
    for src, col in (("li2019_supp1", "medium_raw_li2019"), ("anno2019_growth_medium", "medium_raw_anno2019")):
        vc = prim[col].fillna("<missing>").value_counts()
        m = norm["source"] == src
        norm.loc[m, "n_lines_primary"] = norm.loc[m, "raw_string"].map(vc).fillna(0).astype(int)
    norm["n_lines_primary"] = norm["n_lines_primary"].astype(int)

    both = prim[prim["medium_sources_agree"].notna()]
    agreement = pd.crosstab(both["medium_base_li2019"], both["medium_base_anno2019"])

    inv, inv_info = audit.metabolite_inventory(cells)
    fmt = audit.metabolomics_format(cells)
    expr, colsum = audit.expression_format(cells)

    # acceptance gates
    matched_solid = cells["has_rnaseq"] & cells["in_model"] & ~cells["is_duplicate_line"] & cells["is_solid_tumour"]
    n_primary = int(cells["in_primary_set"].sum())
    coverage = float(cells.loc[matched_solid, "medium_assigned"].mean())
    gates = {
        "n_matched_solid_lines_before_medium": int(matched_solid.sum()),
        "n_primary_lines": n_primary,
        "gate_min_primary_lines": MIN_PRIMARY_LINES,
        "gate_primary_lines_pass": n_primary >= MIN_PRIMARY_LINES,
        "medium_coverage": coverage,
        "gate_min_medium_coverage": MIN_MEDIUM_COVERAGE,
        "gate_medium_coverage_pass": coverage >= MIN_MEDIUM_COVERAGE,
        "medium_coverage_from_li2019_only": float(
            (cells.loc[matched_solid, "medium_source"] == "li2019_supp1").mean()),
        "medium_coverage_from_anno2019_only": float(
            (~cells.loc[matched_solid, "medium_base_anno2019"].isin(medium.UNASSIGNED)).mean()),
        "medium_agreement_n_both": int(len(both)),
        "medium_agreement_frac": float(both["medium_sources_agree"].astype(bool).mean()) if len(both) else None,
    }
    summary = {
        "model_release": res["model_release"],
        "gates": gates,
        "metabolomics_format": fmt,
        "metabolite_inventory": {
            **inv_info,
            "assay_class_counts": inv["assay_class"].value_counts().sort_index().to_dict(),
            "polar_ion_mode_counts": inv.loc[inv["assay_class"].str.startswith("polar"), "ion_mode"]
            .value_counts().sort_index().to_dict(),
            "n_names_with_slash": int(inv["name_has_slash"].sum()),
            "n_polar_carnitines": int((inv["is_acylcarnitine_or_carnitine"] & inv["assay_class"].str.startswith("polar")).sum()),
            "polar_raw_missing_frac_primary_median": float(
                inv.loc[inv["assay_class"].str.startswith("polar"), "raw_missing_frac_primary"].median()),
            "missingness_rule_max_frac": audit.MAX_PREIMPUTATION_MISSING,
            "polar_n_pass_missingness_rule": int(
                inv.loc[inv["assay_class"].str.startswith("polar"), "missingness_rule_pass"].sum()),
            "polar_n_raw_missing_primary_gt_20pct": int(
                (inv.loc[inv["assay_class"].str.startswith("polar"), "raw_missing_frac_primary"] > 0.20).sum()),
        },
        "expression_format": expr,
        "primary_set": {
            "lineage_counts": prim["lineage"].value_counts().sort_index().to_dict(),
            "lineage_collapsed_counts": prim["lineage_collapsed"].value_counts().sort_index().to_dict(),
            "medium_base_counts": prim["medium_base"].value_counts().sort_index().to_dict(),
            "medium_category_counts": prim["medium_category"].value_counts().sort_index().to_dict(),
            "growth_pattern_counts": prim["growth_pattern"].fillna("NA").value_counts().sort_index().to_dict(),
            "n_adherent": int((prim["growth_pattern"] == "Adherent").sum()),
            "n_problematic": int(prim["is_problematic"].sum()),
            "n_sensitivity_excl_problematic": int((~prim["is_problematic"]).sum()),
            "n_in_multi_line_donor_groups": int((prim["donor_group_size_primary"] > 1).sum()),
            "n_multi_line_donor_groups": int(prim.loc[prim["donor_group_size_primary"] > 1, "donor_group_id"].nunique()),
        },
        "issue_counts": issues["issue"].value_counts().sort_index().to_dict(),
    }

    out = {"cells": cells, "issues": issues, "attrition": attrition, "medium_normalization": norm,
           "medium_agreement": agreement, "metabolite_inventory": inv, "summary": _jsonable(summary),
           "rnaseq_colsum": colsum}
    if write:
        _write_csv(cells[CELL_COLUMNS], "cell_lines_phase0.csv")
        _write_csv(norm, "medium_normalization.csv")
        agreement.to_csv(PROCESSED / "medium_source_agreement.csv", lineterminator="\n")
        _write_csv(inv, "metabolite_inventory.csv")
        _write_csv(issues, "id_issues.csv")
        _write_csv(attrition, "attrition.csv")
        (PROCESSED / "audit_summary.json").write_text(json.dumps(out["summary"], indent=2, sort_keys=True) + "\n")
        rows = []
        for f in ["cell_lines_phase0.csv", "medium_normalization.csv", "medium_source_agreement.csv",
                  "metabolite_inventory.csv", "id_issues.csv", "attrition.csv", "audit_summary.json"]:
            sha, _, n = hash_file(PROCESSED / f)
            rows.append({"file": f"data/processed/{f}", "sha256": sha, "bytes": n})
        pd.DataFrame(rows).to_csv(PROCESSED / "FREEZE.tsv", sep="\t", index=False, lineterminator="\n")
    return out


if __name__ == "__main__":
    o = run()
    print(json.dumps(o["summary"]["gates"], indent=2))
    print(f"wrote outputs to {PROCESSED.relative_to(REPO)}")
