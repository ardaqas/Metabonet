"""Metabolite inventory and data-format audits (no modelling, no metabolite->gene mapping)."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from . import load

# Pre-specified target rule (expectations.md §1a, decision D18): a metabolite is eligible only if
# <= 20 % of its values in the primary set were missing before the data generators' MICE imputation.
MAX_PREIMPUTATION_MISSING = 0.20

# Data-column name -> name used in Li et al. '1-LC-MS parameters (polar)' (explicit, hand-checked)
POLAR_ALIASES = {
    "AMP": ["adenosine monophosphate (AMP)"],
    "CMP": ["cytidine monophosphate (CMP)"],
    "GMP": ["guanosine monophosphate (GMP)"],
    "UMP": ["uridine 5'-monophosphate (UMP)"],
    "NADP": ["nicotinamide adenine dinucleotide phosphate (NADP)"],
    "PEP": ["phosphoenolpyruvate (PEP)"],
    "F1P/F6P/G1P/G6P": ["hexose monophosphate (F1P/F6P/G1P/G6P)"],
    # one data feature, two transitions listed in the parameter sheet
    "DHAP/glyceraldehyde 3P": ["dihydroxyacetone phosphate (DHAP)", "glyceraldehyde-3-phosphate"],
}


def _li_row_map(li_index: pd.Index, ccle_ids: pd.Index) -> dict[str, str]:
    ids = set(ccle_ids)
    by_stripped: dict[str, list[str]] = {}
    for c in ccle_ids:
        by_stripped.setdefault(c.split("_", 1)[0], []).append(c)
    out = {}
    for n in li_index:
        if n in ids:
            out[n] = n
        else:
            hits = by_stripped.get(str(n).split("_", 1)[0], [])
            if len(hits) == 1:
                out[n] = hits[0]
    return out


def metabolite_inventory(cells: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    metab = load.load_metabolomics().set_index("CCLE_ID")
    X = metab.drop(columns=["DepMap_ID"])
    polar = load.load_li("lcms_polar")
    lipid = load.load_li("lcms_lipid")
    cv = load.load_li("cv").set_index("Name of metabolites")["CV"]
    lin = load.load_li("lineage_effects").set_index("Name of metabolites")["Lineage effects"]
    raw = load.load_li("raw")
    raw.index = raw.index.map(lambda n: _li_row_map(pd.Index([n]), X.index).get(n, n))

    polar = polar.assign(Name=polar["Name"].astype(str).str.strip()).set_index("Name")
    lipid_names = set(lipid["Lipid"].astype(str).str.strip())
    primary = cells.loc[cells["in_primary_set"], "CCLE_ID"]

    rows = []
    for m in X.columns:
        names = POLAR_ALIASES.get(m, [m])
        in_polar = [n for n in names if n in polar.index]
        if in_polar:
            dps = polar.loc[in_polar, "Declustering potential"].astype(float)
            modes = sorted({"positive" if v > 0 else "negative" for v in dps})
            assay, mode = "polar (HILIC)", "/".join(modes)
            rt = "; ".join(str(polar.at[n, "Retention time"]) for n in in_polar)
            sheet_name = "; ".join(in_polar)
        elif m in lipid_names:
            assay, mode, rt, sheet_name = "lipid (RPLC)", "", "", m
        else:
            assay, mode, rt, sheet_name = "UNCLASSIFIED", "", "", ""
        r = raw[m] if m in raw.columns else pd.Series(dtype=float)
        x_all, x_pri = X[m], X.loc[X.index.intersection(primary), m]
        r_pri = r.reindex(primary)
        rows.append({
            "metabolite": m,
            "assay_class": assay,
            "ion_mode": mode,
            "lcms_sheet_name": sheet_name,
            "retention_time_min": rt,
            "name_has_slash": "/" in m,
            "is_acylcarnitine_or_carnitine": bool(re.search(r"carnitine", m, re.I)),
            "li2019_cv": cv.get(m, np.nan),
            "li2019_lineage_r2": lin.get(m, np.nan),
            "raw_missing_frac_all": float(r.isna().mean()) if len(r) else np.nan,
            "raw_missing_frac_primary": float(r_pri.isna().mean()) if len(r) else np.nan,
            "n_missing_released_csv": int(x_all.isna().sum()),
            "mean_primary": float(x_pri.mean()),
            "sd_primary": float(x_pri.std()),
            "min_primary": float(x_pri.min()),
            "max_primary": float(x_pri.max()),
            "n_unique_primary": int(x_pri.nunique()),
        })
    inv = pd.DataFrame(rows)
    inv["missingness_rule_pass"] = inv["raw_missing_frac_primary"] <= MAX_PREIMPUTATION_MISSING
    used = {n for v in POLAR_ALIASES.values() for n in v} | set(X.columns)
    info = {
        "polar_sheet_names_unused": sorted(set(polar.index) - used),
        "n_polar_transitions_in_sheet": int(len(polar)),
        "n_lipid_rows_in_sheet": int(len(lipid)),
    }
    return inv, info


def metabolomics_format(cells: pd.DataFrame) -> dict:
    """Compare released CSV to Li clean data; infer transformation from raw peak areas."""
    metab = load.load_metabolomics().set_index("CCLE_ID").drop(columns=["DepMap_ID"])
    clean = load.load_li("clean")
    raw = load.load_li("raw")
    cmap = _li_row_map(clean.index, metab.index)
    rmap = _li_row_map(raw.index, metab.index)
    clean.index = clean.index.map(lambda n: cmap.get(n, n))
    raw.index = raw.index.map(lambda n: rmap.get(n, n))
    common_r = metab.index.intersection(clean.index)
    common_c = metab.columns.intersection(clean.columns)
    diff = (metab.loc[common_r, common_c] - clean.loc[common_r, common_c].astype(float)).abs()

    rr = metab.index.intersection(raw.index)
    cc = metab.columns.intersection(raw.columns)
    R = raw.loc[rr, cc].apply(pd.to_numeric, errors="coerce")
    C = metab.loc[rr, cc]
    corr, slope, lin_corr = [], [], []
    for m in cc:
        ok = R[m].notna() & (R[m] > 0)
        if ok.sum() > 10:
            lr = np.log10(R.loc[ok, m].astype(float))
            corr.append(np.corrcoef(lr, C.loc[ok, m])[0, 1])
            slope.append(np.polyfit(lr, C.loc[ok, m], 1)[0])
            lin_corr.append(np.corrcoef(R.loc[ok, m].astype(float), C.loc[ok, m])[0, 1])
    imputed_cells = int(R.isna().sum().sum())
    vals = metab.values.ravel()
    vals = vals[~np.isnan(vals)]
    return {
        "csv_shape": metab.shape,
        "clean_sheet_shape": clean.shape,
        "rows_matched_csv_vs_clean": len(common_r),
        "cols_matched_csv_vs_clean": len(common_c),
        "max_abs_diff_csv_vs_clean": float(np.nanmax(diff.values)),
        "raw_rows_matched": len(rr),
        "raw_missing_cells": imputed_cells,
        "raw_missing_frac": imputed_cells / R.size,
        "raw_missing_filled_in_csv": int(C[R.isna()].notna().sum().sum()),
        "median_corr_clean_vs_log10raw": float(np.median(corr)),
        "min_corr_clean_vs_log10raw": float(np.min(corr)),
        "median_corr_clean_vs_linear_raw": float(np.median(lin_corr)),
        "median_slope_clean_on_log10raw": float(np.median(slope)),
        "iqr_slope": (float(np.percentile(slope, 25)), float(np.percentile(slope, 75))),
        "value_min": float(vals.min()), "value_p01": float(np.percentile(vals, 1)),
        "value_median": float(np.median(vals)), "value_p99": float(np.percentile(vals, 99)),
        "value_max": float(vals.max()),
        "n_metabolites_corr_computed": len(corr),
    }


def expression_format(cells: pd.DataFrame) -> tuple[dict, pd.Series]:
    """Audit the CCLE 2019 RSEM TPM matrix. Returns (summary, per-line TPM column sums)."""
    tpm, tx = load.load_rnaseq()
    ids = tpm.index.to_series()
    stripped = ids.str.replace(r"\.\d+$", "", regex=True)
    is_r = ids.str.startswith("ENSGR")
    r_counterpart = stripped[is_r].str.replace("ENSGR", "ENSG0", regex=False).isin(set(stripped[~is_r]))
    primary = cells.loc[cells["in_primary_set"], "CCLE_ID"]
    P = tpm[primary]
    logP = np.log2(P + 1.0)
    colsum = tpm.sum(axis=0)
    mean_log = logP.mean(axis=1)
    summary = {
        "n_genes": int(tpm.shape[0]),
        "n_lines": int(tpm.shape[1]),
        "n_ensg_versioned_ids": int(ids.str.match(r"^ENSG\d{11}\.\d+$").sum()),
        "n_ensgr_par_y_ids": int(ids.str.match(r"^ENSGR\d{10}\.\d+$").sum()),
        "n_other_ids": int((~ids.str.match(r"^ENSGR?\d{10,11}\.\d+$")).sum()),
        "n_duplicate_ids_after_version_strip": int(stripped.duplicated().sum()),
        "n_nan": int(np.isnan(tpm.values).sum()),
        "n_negative": int((tpm.values < 0).sum()),
        "max_value": float(np.nanmax(tpm.values)),
        "colsum_min": float(colsum.min()), "colsum_median": float(colsum.median()), "colsum_max": float(colsum.max()),
        "frac_zero": float((tpm.values == 0).mean()),
        "n_primary_lines": int(P.shape[1]),
        "n_genes_mean_log2tpm1_ge1_primary": int((mean_log >= 1).sum()),
        "n_genes_all_zero_primary": int((P.sum(axis=1) == 0).sum()),
        "duplicate_line_columns": int(pd.Index(tpm.columns).duplicated().sum()),
        "n_ensgr_with_ensg_counterpart": int(r_counterpart.sum()),
        "n_ensgr_all_zero_all_lines": int((tpm[is_r.values].sum(axis=1) == 0).sum()),
        "n_genes_after_dropping_ensgr": int((~is_r).sum()),
    }
    return summary, colsum
