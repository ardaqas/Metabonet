"""Typed loaders for Phase 0 raw files. Raw files are read as data only."""
from __future__ import annotations

import re
import warnings
from functools import lru_cache

import numpy as np
import openpyxl
import pandas as pd

from .sources import MODEL_RELEASE_FILE, SOURCES

LI_SHEETS = {
    "annotations": "1-cell line annotations",
    "clean": "1-clean data",
    "raw": "1-raw data",
    "vendors": "1-media vendors",
    "cv": "1-CV",
    "lineage_effects": "1-lineage effects",
    "lcms_polar": "1-LC-MS parameters (polar)",
    "lcms_lipid": "1-LC-MS parameters (lipid)",
}

MODEL_REQUIRED = ["ModelID", "OncotreeLineage", "RRID", "PatientID", "GrowthPattern"]


def _path(key: str):
    p = SOURCES[key].path
    if not p.exists():
        raise FileNotFoundError(f"{p} not found; run `make data` (see docs/data_sources.md)")
    return p


# --------------------------------------------------------------------------- CCLE 2019
def load_metabolomics() -> pd.DataFrame:
    """928 lines x (CCLE_ID, DepMap_ID, 225 metabolites). Returned as-is (no transformation)."""
    return pd.read_csv(_path("metabolomics"))


def load_annotations() -> pd.DataFrame:
    return pd.read_csv(_path("annotations"), sep="\t", dtype=str)


def load_rnaseq_header() -> list[str]:
    """Column names of the RSEM TPM file without reading the matrix."""
    return list(pd.read_csv(_path("rnaseq_tpm"), sep="\t", nrows=0).columns)


def load_rnaseq() -> tuple[pd.DataFrame, pd.Series]:
    """Return (TPM matrix genes x CCLE_IDs [float32], transcript_ids Series)."""
    df = pd.read_csv(_path("rnaseq_tpm"), sep="\t", index_col=0, low_memory=False)
    tx = df.pop("transcript_ids")
    return df.astype(np.float32), tx


def drop_par_y(tpm: pd.DataFrame) -> pd.DataFrame:
    """Drop GENCODE v19 chrY pseudoautosomal copies (ENSGR... IDs; decision D22).

    In the CCLE 2019 RSEM matrix all 47 ENSGR rows are zero in every line (reads are
    assigned to the chrX copy, which keeps its ENSG... ID), so nothing is lost.
    """
    return tpm.loc[~tpm.index.astype(str).str.startswith("ENSGR")]


# --------------------------------------------------------------------------- Li et al. 2019
@lru_cache(maxsize=1)
def _li_book() -> pd.ExcelFile:
    warnings.filterwarnings("ignore", message="Unknown extension", module="openpyxl")
    return pd.ExcelFile(_path("li2019_supp1"), engine="openpyxl")


def load_li(sheet: str) -> pd.DataFrame:
    """Load one Li et al. Supp. Table 1 sheet by short key (see LI_SHEETS)."""
    name = LI_SHEETS[sheet]
    if sheet in ("lcms_polar", "lcms_lipid"):
        return _read_until_blank(name)
    if sheet in ("clean", "raw"):
        df = pd.read_excel(_li_book(), sheet_name=name, index_col=0, na_values=["NA"])
        df.index.name = "CCLE_ID"
        return df
    return pd.read_excel(_li_book(), sheet_name=name)


def _read_until_blank(sheet_name: str) -> pd.DataFrame:
    """The LC-MS parameter sheets declare ~1M (empty) rows; stop at the first blank row."""
    warnings.filterwarnings("ignore", message="Unknown extension", module="openpyxl")
    wb = openpyxl.load_workbook(_path("li2019_supp1"), read_only=True)
    ws = wb[sheet_name]
    rows = []
    for row in ws.iter_rows(values_only=True):
        if all(v is None for v in row):
            break
        rows.append(row)
    wb.close()
    header = [str(h).strip() if h is not None else f"col{i}" for i, h in enumerate(rows[0])]
    body = rows[1:]
    # polar sheet has a units row: (m/z), (DP), ...
    if body and body[0][0] is None:
        body = body[1:]
    df = pd.DataFrame(body, columns=header)
    return df.loc[:, [c for c in df.columns if not (c.startswith("col") and df[c].isna().all())]]


# --------------------------------------------------------------------------- DepMap
def release_evidence() -> dict:
    """Evidence for which DepMap release Model.csv belongs to.

    * ``readme_release``: header line of data/raw/depmap/README.txt (e.g. 'DepMap Public 26Q3')
    * ``release_txt``: first line of data/raw/depmap/RELEASE.txt (written at download time)
    * ``schema_match``: Model.csv columns == columns documented in README's '### Model.csv' section (same order)
    """
    readme, model = SOURCES["model_readme"].path, SOURCES["model"].path
    ev = {"readme_release": None, "release_txt": None, "schema_match": None,
          "n_documented_columns": None, "n_file_columns": None}
    if MODEL_RELEASE_FILE.exists():
        ev["release_txt"] = MODEL_RELEASE_FILE.read_text().strip().splitlines()[0].strip()
    if readme.exists():
        text = readme.read_text(errors="replace")
        m = re.match(r"\s*#\s*(DepMap Public \S+)", text)
        ev["readme_release"] = m.group(1) if m else None
        if "### Model.csv" in text and model.exists():
            sec = text[text.index("### Model.csv"):]
            nxt = sec.find("\n### ", 5)
            sec = sec if nxt < 0 else sec[:nxt]
            doc = re.findall(r"^- (\w+):", sec, re.M)
            cols = list(pd.read_csv(model, nrows=0).columns)
            ev.update(schema_match=(doc == cols), n_documented_columns=len(doc), n_file_columns=len(cols))
    return ev


def model_release() -> str:
    ev = release_evidence()
    if ev["readme_release"]:
        schema = "Model.csv schema matches README" if ev["schema_match"] else "Model.csv schema does NOT match README"
        return f"{ev['readme_release']} (from depmap/README.txt header; {schema})"
    if ev["release_txt"]:
        return f"{ev['release_txt']} (from RELEASE.txt only)"
    return "UNKNOWN (no README.txt or RELEASE.txt in data/raw/depmap/)"


def load_model() -> tuple[pd.DataFrame, list[str]]:
    """Return (Model.csv, list of notes about schema fallbacks)."""
    m = pd.read_csv(_path("model"), dtype=str)
    notes = []
    missing = [c for c in MODEL_REQUIRED if c not in m.columns]
    if missing:
        raise ValueError(f"Model.csv lacks required columns: {missing}")
    if "CCLEName" not in m.columns:
        m["CCLEName"] = np.nan
        notes.append("Model.csv has no CCLEName column; name cross-check uses StrippedCellLineName only")
    if m["ModelID"].duplicated().any():
        raise ValueError("duplicated ModelID in Model.csv")
    return m, notes
