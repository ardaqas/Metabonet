"""Human-GEM (pinned v2.0.0) loading: SBML via COBRApy (parsing only - no FBA, no solver use)
plus the annotation tables metabolites.tsv / genes.tsv / reactions.tsv.
"""
from __future__ import annotations

import pickle
import re
from functools import lru_cache

import pandas as pd

from ..phase0.download import hash_file
from ..phase0.sources import REPO, SOURCES

INTERIM = REPO / "data" / "interim"

# Reaction exclusion rule lives in genesets.EXCLUDED_SUBSYSTEMS (D23, D33)


def _p(name: str):
    return SOURCES[f"humangem_{name}"].path


@lru_cache(maxsize=1)
def load_model():
    """cobra.Model for Human-GEM.xml (cached as pickle keyed on the SBML sha256)."""
    import cobra

    sbml = _p("human_gem_xml")
    sha = hash_file(sbml)[0][:16]
    INTERIM.mkdir(parents=True, exist_ok=True)
    cache = INTERIM / f"Human-GEM_{sha}.pkl"
    if cache.exists():
        return pickle.loads(cache.read_bytes())
    model = cobra.io.read_sbml_model(str(sbml))
    cache.write_bytes(pickle.dumps(model))
    return model


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(_p(f"{name}_tsv"), sep="\t", dtype=str).fillna("")


def split_ids(s: str) -> list[str]:
    return [x.strip() for x in re.split(r"[;,]", s or "") if x.strip()]


def norm_hmdb(h: str) -> str:
    """HMDB00123 -> HMDB0000123 (HMDB v4 7-digit form)."""
    m = re.fullmatch(r"HMDB(\d+)", h.strip())
    return f"HMDB{int(m.group(1)):07d}" if m else h.strip()


def norm_chebi(c: str) -> str:
    c = c.strip()
    return c if c.startswith("CHEBI:") else (f"CHEBI:{c}" if c.isdigit() else c)


def metabolite_table() -> pd.DataFrame:
    """One row per compartment-free metabolite (metsNoComp) with names and external IDs."""
    model = load_model()
    tsv = read_tsv("metabolites")
    names = {m.id: m.name for m in model.metabolites}
    formulas = {m.id: m.formula for m in model.metabolites}
    tsv["name"] = tsv["mets"].map(names)
    tsv["formula"] = tsv["mets"].map(formulas)
    tsv["comp"] = tsv["mets"].str[-1]
    rows = []
    for nc, g in tsv.groupby("metsNoComp"):
        rows.append({
            "metsNoComp": nc,
            "name": g["name"].iloc[0],
            "names_all": "; ".join(sorted(set(g["name"]))),
            "formula": g["formula"].iloc[0],
            "compartments": "".join(sorted(g["comp"])),
            "n_copies": len(g),
            "kegg": sorted({k for s in g["metKEGGID"] for k in split_ids(s)}),
            "hmdb": sorted({norm_hmdb(k) for s in g["metHMDBID"] for k in split_ids(s)}),
            "chebi": sorted({norm_chebi(k) for s in g["metChEBIID"] for k in split_ids(s)}),
            "bigg": sorted({k for s in g["metBiGGID"] for k in split_ids(s)}),
            "metanetx": sorted({k for s in g["metMetaNetXID"] for k in split_ids(s)}),
        })
    return pd.DataFrame(rows)


def gene_table() -> pd.DataFrame:
    return read_tsv("genes")
