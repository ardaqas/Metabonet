"""Cell-line crosswalk and the Phase 0 matched cell-line table (research_plan.md §4)."""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import cellosaurus, load, medium

HAEM_LINEAGES = {"Lymphoid", "Myeloid"}
MIN_LINEAGE_N = 15


@dataclass
class Issues:
    rows: list[dict] = field(default_factory=list)

    def add(self, ccle_id, issue, detail, action, depmap_id=None):
        self.rows.append({"CCLE_ID": ccle_id, "DepMap_ID": depmap_id, "issue": issue,
                          "detail": detail, "action": action})

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows, columns=["CCLE_ID", "DepMap_ID", "issue", "detail", "action"])


def _stripped(ccle_id: str) -> str:
    return ccle_id.split("_", 1)[0]


def _norm_name(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def match_li_names(li_names: pd.Index, ccle_ids: pd.Index, issues: Issues, sheet: str) -> dict[str, str]:
    """Map Li et al. row names -> CCLE_ID (exact, else unique stripped-name match)."""
    ids = set(ccle_ids)
    by_stripped = pd.Series(list(ccle_ids), index=[_stripped(c) for c in ccle_ids])
    out = {}
    for n in li_names:
        if n in ids:
            out[n] = n
            continue
        hits = by_stripped.get(_stripped(n))
        hits = [] if hits is None else ([hits] if isinstance(hits, str) else list(hits))
        if len(hits) == 1:
            out[n] = hits[0]
            issues.add(hits[0], f"li2019_name_differs[{sheet}]",
                       f"Li et al. name '{n}' vs CCLE 2019 '{hits[0]}'", "matched by cell line name (unique)")
        else:
            issues.add(None, f"li2019_name_unmatched[{sheet}]", f"'{n}' ({len(hits)} candidates)", "dropped")
    return out


def build(cello_paths=None) -> dict[str, pd.DataFrame]:
    """Build the matched cell-line table. Returns dict of DataFrames."""
    issues = Issues()
    attr: list[tuple[str, int]] = []

    metab = load.load_metabolomics()
    anno = load.load_annotations()
    model, model_notes = load.load_model()
    rna_cols = load.load_rnaseq_header()[2:]
    li_ann = load.load_li("annotations")

    for note in model_notes:
        issues.add(None, "model_schema", note, "fallback used")

    t = metab[["CCLE_ID", "DepMap_ID"]].rename(columns={"DepMap_ID": "DepMap_ID_metab"}).copy()
    attr.append(("CCLE 2019 metabolomics lines", len(t)))
    for col in ("CCLE_ID",):
        d = t[col][t[col].duplicated()]
        for v in d:
            issues.add(v, "duplicate_ccle_id_metabolomics", "", "investigate")

    # ---- DepMap ID: metabolomics file, cross-checked with 2019 annotation
    a = anno.set_index("CCLE_ID")
    t["DepMap_ID_anno2019"] = t["CCLE_ID"].map(a["depMapID"])
    t["in_anno2019"] = t["CCLE_ID"].isin(a.index)
    for r in t.itertuples():
        if not r.in_anno2019:
            issues.add(r.CCLE_ID, "missing_in_anno2019", "", "medium fallback unavailable", r.DepMap_ID_metab)
        if pd.isna(r.DepMap_ID_metab) and pd.notna(r.DepMap_ID_anno2019):
            issues.add(r.CCLE_ID, "depmap_id_missing_in_metabolomics",
                       f"filled from 2019 annotation: {r.DepMap_ID_anno2019}", "filled", r.DepMap_ID_anno2019)
        elif pd.notna(r.DepMap_ID_metab) and pd.notna(r.DepMap_ID_anno2019) and r.DepMap_ID_metab != r.DepMap_ID_anno2019:
            issues.add(r.CCLE_ID, "depmap_id_conflict_metab_vs_anno",
                       f"{r.DepMap_ID_metab} vs {r.DepMap_ID_anno2019}", "metabolomics ID kept", r.DepMap_ID_metab)
    t["DepMap_ID"] = t["DepMap_ID_metab"].fillna(t["DepMap_ID_anno2019"])
    for v in t["DepMap_ID"][t["DepMap_ID"].duplicated() & t["DepMap_ID"].notna()]:
        issues.add(None, "duplicate_depmap_id", v, "investigate", v)
    attr.append(("... with a DepMap ID", int(t["DepMap_ID"].notna().sum())))

    # ---- Model.csv (current DepMap): ModelID, else ModelIDAlias
    mi = model.set_index("ModelID")
    alias = {}
    if "ModelIDAlias" in model.columns:
        for mid, al in zip(model["ModelID"], model["ModelIDAlias"]):
            if isinstance(al, str):
                for x in re.split(r"[;,\s]+", al.strip()):
                    if x:
                        alias[x] = mid
    t["ModelID"] = np.where(t["DepMap_ID"].isin(mi.index), t["DepMap_ID"], None)
    for i, r in t[t["ModelID"].isna() & t["DepMap_ID"].notna()].iterrows():
        if r.DepMap_ID in alias:
            t.at[i, "ModelID"] = alias[r.DepMap_ID]
            issues.add(r.CCLE_ID, "depmap_id_retired_alias",
                       f"{r.DepMap_ID} is listed as ModelIDAlias of {alias[r.DepMap_ID]} in current Model.csv",
                       "mapped to current ModelID", r.DepMap_ID)
        else:
            issues.add(r.CCLE_ID, "not_in_model_csv", r.DepMap_ID, "no lineage; excluded", r.DepMap_ID)
    t["in_model"] = t["ModelID"].notna()
    attr.append(("... found in current Model.csv (ModelID or alias)", int(t["in_model"].sum())))

    # duplicate cell line: DepMap merged this ID into a ModelID that is itself a metabolomics row
    t["is_duplicate_line"] = False
    canon = set(t.loc[t["DepMap_ID"] == t["ModelID"], "ModelID"])
    for i, r in t[(t["ModelID"] != t["DepMap_ID"]) & t["ModelID"].notna()].iterrows():
        if r.ModelID in canon:
            t.at[i, "is_duplicate_line"] = True
            other = t.loc[(t["DepMap_ID"] == r.ModelID), "CCLE_ID"].iloc[0]
            issues.add(r.CCLE_ID, "duplicate_cell_line",
                       f"same current DepMap model ({r.ModelID}) as metabolomics line {other}",
                       f"excluded from primary set; canonical line {other} kept", r.DepMap_ID)

    mcols = ["CCLEName", "StrippedCellLineName", "CellLineName", "RRID", "PatientID", "OncotreeLineage",
             "OncotreePrimaryDisease", "OncotreeSubtype", "DepmapModelType", "GrowthPattern", "FormulationID"]
    mcols = [c for c in mcols if c in mi.columns]
    t = t.join(mi[mcols], on="ModelID")
    for r in t[t["in_model"]].itertuples():
        if isinstance(r.CCLEName, str) and r.CCLEName != r.CCLE_ID:
            same_stripped = _stripped(r.CCLEName) == _stripped(r.CCLE_ID)
            issues.add(r.CCLE_ID, "ccle_name_changed_in_model_csv",
                       f"Model.csv CCLEName='{r.CCLEName}' ({'tissue suffix only' if same_stripped else 'different cell line name'})",
                       "kept (DepMap ID consistent); lineage from Model.csv", r.DepMap_ID)
        elif not isinstance(r.CCLEName, str) and _norm_name(r.StrippedCellLineName) != _norm_name(_stripped(r.CCLE_ID)):
            issues.add(r.CCLE_ID, "stripped_name_differs_in_model_csv", str(r.StrippedCellLineName), "kept", r.DepMap_ID)

    # ---- RNA-seq (CCLE 2019 RSEM TPM; columns are CCLE_IDs)
    t["has_rnaseq"] = t["CCLE_ID"].isin(set(rna_cols))
    for r in t[~t["has_rnaseq"]].itertuples():
        issues.add(r.CCLE_ID, "no_ccle2019_rnaseq", "CCLE_ID absent from RSEM TPM columns (no alias match)",
                   "excluded", r.DepMap_ID)
    t["has_metabolomics"] = True
    attr.append(("... with CCLE 2019 RNA-seq", int((t["in_model"] & t["has_rnaseq"]).sum())))

    # ---- Lineage, haematopoietic and non-cancer flags
    t["lineage"] = t["OncotreeLineage"]
    t["site_primary_2019"] = t["CCLE_ID"].map(a["Site_Primary"])
    t["is_haematopoietic"] = t["lineage"].isin(HAEM_LINEAGES)
    haem2019 = t["site_primary_2019"].eq("haematopoietic_and_lymphoid_tissue") | t["CCLE_ID"].str.endswith(
        "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")
    for r in t[(t["is_haematopoietic"] != haem2019) & t["lineage"].notna()].itertuples():
        issues.add(r.CCLE_ID, "haematopoietic_flag_discordant",
                   f"Model.csv lineage={r.lineage}; 2019 Site_Primary={r.site_primary_2019}",
                   "Model.csv lineage used", r.DepMap_ID)
    t["is_non_cancerous"] = t["OncotreePrimaryDisease"].eq("Non-Cancerous") | t["lineage"].eq("Fibroblast")
    for r in t[t["is_non_cancerous"]].itertuples():
        issues.add(r.CCLE_ID, "non_cancerous_model",
                   f"lineage={r.lineage}; disease={r.OncotreePrimaryDisease}; type={r.DepmapModelType}",
                   "excluded from solid-tumour set", r.DepMap_ID)
    t["growth_pattern"] = t["GrowthPattern"]
    t["is_solid_tumour"] = t["lineage"].notna() & ~t["is_haematopoietic"] & ~t["is_non_cancerous"]

    # ---- Medium (Li et al. preferred; 2019 Growth.Medium fallback)
    li_map = match_li_names(li_ann["Name with tissue origins"], pd.Index(t["CCLE_ID"]), issues, "annotations")
    li_med = li_ann.assign(CCLE_ID=li_ann["Name with tissue origins"].map(li_map)).dropna(subset=["CCLE_ID"])
    li_med = li_med.set_index("CCLE_ID")
    ti = t.set_index("CCLE_ID")
    med = medium.assign(li_med["Cell culture media"].reindex(ti.index), a["Growth.Medium"].reindex(ti.index))
    t = t.join(med, on="CCLE_ID")
    t["li2019_classification"] = t["CCLE_ID"].map(li_med["Classifications"])
    t["medium_assigned"] = ~t["medium_base"].isin(medium.UNASSIGNED)

    # ---- Cellosaurus: problematic lines and donor relationships
    t = _add_cellosaurus(t, issues, cello_paths)

    # ---- Primary set
    t["in_primary_set"] = (t["has_metabolomics"] & t["has_rnaseq"] & t["in_model"] & ~t["is_duplicate_line"]
                           & t["is_solid_tumour"] & t["medium_assigned"])
    step = t["has_rnaseq"] & t["in_model"] & t["lineage"].notna()
    attr.append(("... with lineage (OncotreeLineage)", int(step.sum())))
    step &= ~t["is_duplicate_line"]
    attr.append(("... duplicate cell lines removed (DepMap-merged IDs)", int(step.sum())))
    step &= ~t["is_non_cancerous"]
    attr.append(("... cancer models (non-cancerous/fibroblast removed)", int(step.sum())))
    step &= ~t["is_haematopoietic"]
    attr.append(("... solid tumour (Lymphoid/Myeloid removed)", int(step.sum())))
    step &= t["medium_assigned"]
    attr.append(("... with medium assigned = PRIMARY SET", int(step.sum())))
    assert step.equals(t["in_primary_set"])

    # ---- Category collapsing on the primary set (< 15 lines -> other)
    t["lineage_collapsed"] = _collapse_lineage(t["lineage"], t["in_primary_set"])
    t["medium_category"], _ = medium.collapse(t["medium_base"], t["in_primary_set"])

    # donor groups among primary set
    t["donor_group_size_primary"] = t["donor_group_id"].map(t.loc[t["in_primary_set"], "donor_group_id"].value_counts()).fillna(0).astype(int)

    attrition = pd.DataFrame(attr, columns=["step", "n_lines"])
    return {"cells": t, "issues": issues.frame(), "attrition": attrition, "model_release": load.model_release()}


def _collapse_lineage(lin: pd.Series, in_set: pd.Series) -> pd.Series:
    counts = lin[in_set].value_counts()
    keep = set(counts[counts >= MIN_LINEAGE_N].index)
    return lin.where(lin.isin(keep) | lin.isna(), "Other")


class _UF:
    def __init__(self, items):
        self.p = {i: i for i in items}

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def _add_cellosaurus(t: pd.DataFrame, issues: Issues, paths) -> pd.DataFrame:
    if paths is None:
        paths = sorted(cellosaurus.CELLOSAURUS_DIR.glob("batch_*.txt"))
    if not paths:
        raise FileNotFoundError("No Cellosaurus cache; run `python -m src.phase0.cellosaurus` (make data)")
    cdf = cellosaurus.load_records(paths)
    idx = cellosaurus.accession_index(cdf)
    t["cellosaurus_ac"] = t["RRID"].map(lambda r: idx.get(r) if isinstance(r, str) else None)
    for r in t[t["RRID"].notna() & t["cellosaurus_ac"].isna()].itertuples():
        issues.add(r.CCLE_ID, "rrid_not_found_in_cellosaurus", r.RRID, "no problematic/donor info", r.DepMap_ID)
    for r in t[(t["cellosaurus_ac"].notna()) & (t["cellosaurus_ac"] != t["RRID"])].itertuples():
        issues.add(r.CCLE_ID, "rrid_secondary_accession", f"{r.RRID} -> {r.cellosaurus_ac}", "primary used", r.DepMap_ID)
    c = cdf.set_index("AC")
    t["cellosaurus_name"] = t["cellosaurus_ac"].map(c["ID"])
    t["is_problematic"] = t["cellosaurus_ac"].map(c["is_problematic"]).fillna(False).astype(bool)
    t["problematic_type"] = t["cellosaurus_ac"].map(c["problematic_type"]).fillna("")
    t["problematic_text"] = t["cellosaurus_ac"].map(c["problematic_text"]).fillna("")
    for r in t[t["is_problematic"]].itertuples():
        issues.add(r.CCLE_ID, "cellosaurus_problematic", r.problematic_text[:300],
                   "flagged; kept (see decision log)", r.DepMap_ID)

    # duplicates: same RRID on two metabolomics rows
    for rrid, g in t[t["cellosaurus_ac"].notna()].groupby("cellosaurus_ac"):
        if len(g) > 1:
            for cid in g["CCLE_ID"]:
                issues.add(cid, "duplicate_rrid", f"{rrid} shared by {list(g['CCLE_ID'])}",
                           "same donor group", None)

    # donor groups: PatientID, same RRID, Cellosaurus OI / HI, 'shown to be X derivative'
    uf = _UF(list(t["CCLE_ID"]))
    reason: dict[str, set] = {cid: set() for cid in t["CCLE_ID"]}

    def link(a, b, why):
        if a != b:
            uf.union(a, b)
            reason[a].add(why)
            reason[b].add(why)

    for key, why in (("PatientID", "same DepMap PatientID"), ("cellosaurus_ac", "same RRID")):
        for _, g in t[t[key].notna()].groupby(key):
            ids = list(g["CCLE_ID"])
            for x in ids[1:]:
                link(ids[0], x, why)
    ac2cid = {}
    for cid, ac in zip(t["CCLE_ID"], t["cellosaurus_ac"]):
        if isinstance(ac, str):
            ac2cid.setdefault(ac, cid)
    name2cid = {}
    for cid, ac in zip(t["CCLE_ID"], t["cellosaurus_ac"]):
        if isinstance(ac, str):
            name2cid.setdefault(_norm_name(c.at[ac, "ID"]), cid)
            for s in c.at[ac, "SY"]:
                name2cid.setdefault(_norm_name(s), cid)
    for cid, ac in zip(t["CCLE_ID"], t["cellosaurus_ac"]):
        if not isinstance(ac, str):
            continue
        for o in c.at[ac, "OI"]:
            if idx.get(o) in ac2cid:
                link(cid, ac2cid[idx[o]], "Cellosaurus OI (same individual)")
        for h in c.at[ac, "HI"]:
            if idx.get(h) in ac2cid:
                link(cid, ac2cid[idx[h]], "Cellosaurus HI (derived from)")
        dn = c.at[ac, "derivative_of_name"]
        if dn and _norm_name(dn) in name2cid:
            link(cid, name2cid[_norm_name(dn)], "Cellosaurus problematic: shown to be derivative")

    roots = {cid: uf.find(cid) for cid in t["CCLE_ID"]}
    sizes = pd.Series(roots).value_counts()
    order = sorted({r for r in roots.values() if sizes[r] > 1})
    gid = {r: f"DG{i + 1:03d}" for i, r in enumerate(order)}
    t["donor_group_id"] = t["CCLE_ID"].map(lambda x: gid.get(roots[x], f"S_{x}"))
    t["donor_group_size"] = t["CCLE_ID"].map(lambda x: int(sizes[roots[x]]))
    t["donor_link_reason"] = t["CCLE_ID"].map(lambda x: "; ".join(sorted(reason[x])))
    for r in t[t["donor_group_size"] > 1].itertuples():
        members = list(t.loc[t["donor_group_id"] == r.donor_group_id, "CCLE_ID"])
        issues.add(r.CCLE_ID, "shared_donor_group", f"{r.donor_group_id}: {members} ({r.donor_link_reason})",
                   "keep in same CV fold", r.DepMap_ID)
    return t
