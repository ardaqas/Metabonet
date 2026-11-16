"""Frozen Phase 2 inputs for the primary set (D40).

Everything here is read from Phase 0/1 frozen tables; nothing is recomputed from raw files except the
log2(TPM+1) matrix (via the Phase 1 loader, same definition as rule 3).

Gene universes (D40):
  U_all   = expressed genes (mean log2(TPM+1) >= 1 over the 735 primary lines) that are protein-coding in
            GENCODE v19, plus every expressed gene that maps to a Human-GEM gene (so each target's own set
            is inside the universe). Used for M2 and as the N-all pool.
  U_met   = expressed genes that map to a Human-GEM gene (2,123). N-met pool.
  G_m     = frozen first-neighbour set of target m, expressed genes only (gene_set_members.csv).
Random pools exclude G_m per target.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from ..phase0 import load
from ..phase0.sources import PROCESSED
from ..phase1 import genemap

P1 = PROCESSED / "phase1"


@dataclass
class Inputs:
    cells: pd.DataFrame            # primary lines, ordered; lineage_collapsed, medium_category, donor_group_id
    Y: pd.DataFrame                # lines x features (log10 metabolomics), targets + diagnostics
    X: np.ndarray                  # lines x genes (log2(TPM+1), float32), genes = universe
    genes: pd.Index                # universe gene IDs (CCLE versioned Ensembl)
    gene_type: pd.Series           # v19 gene_type per universe gene ('' if absent from GTF)
    targets: list[str]             # frozen target features (74)
    diagnostics: list[str]         # pipeline diagnostics (MNA)
    sets_primary: dict[str, list[str]]
    sets_enzymatic: dict[str, list[str]]
    met_pool: list[str]            # U_met
    classes: pd.Series             # target class per feature

    @property
    def gene_index(self) -> dict[str, int]:
        return {g: i for i, g in enumerate(self.genes)}

    def idx(self, genes: list[str]) -> np.ndarray:
        gi = self.gene_index
        return np.array([gi[g] for g in genes], dtype=int)


POPULATIONS = {                      # D46: analysis populations; categories re-collapsed (< 15 -> Other) per population
    "primary": lambda c: c["in_primary_set"],
    "all_lines": lambda c: c["has_rnaseq"] & c["in_model"] & ~c["is_duplicate_line"] & ~c["is_non_cancerous"] & c["medium_assigned"],
    "no_problematic": lambda c: c["in_primary_set"] & ~c["is_problematic"],
}
MIN_CATEGORY_N = 15


def _collapse(s: pd.Series, min_n: int = MIN_CATEGORY_N, other: str = "Other") -> pd.Series:
    vc = s.value_counts()
    return s.where(s.isin(vc[vc >= min_n].index), other)


def population_cells(population: str) -> pd.DataFrame:
    cells = pd.read_csv(PROCESSED / "cell_lines_phase0.csv")
    cells = cells[POPULATIONS[population](cells)].sort_values("CCLE_ID").reset_index(drop=True)
    if population != "primary":
        cells["lineage_collapsed"] = _collapse(cells["lineage"])
        base = cells["medium_base"].where(cells["medium_base"] != "mixture", "other")   # D13: mixtures -> other
        cells["medium_category"] = _collapse(base, other="other")
    return cells


@lru_cache(maxsize=1)
def _universe() -> tuple[list[str], list[str], pd.Series]:
    """Gene universe frozen on the primary set (D40); reused unchanged by every population."""
    cells = population_cells("primary")
    lines = cells["CCLE_ID"].tolist()
    eg, _ = genemap.expression_genes(lines)
    gtf = genemap.read_gtf_genes().set_index("gene_id")["gene_type"]
    gmap = pd.read_csv(P1 / "gene_mapping.csv")
    mapped = set(gmap.loc[gmap["mapped"], "ccle_gene_id"])
    expressed = set(eg.loc[eg["expressed"], "gene_id"])
    gtype = eg["gene_id"].map(gtf).fillna("")
    coding = set(eg.loc[gtype.eq("protein_coding").values, "gene_id"])
    universe = sorted((expressed & coding) | (expressed & mapped))
    met_pool = sorted(expressed & mapped)
    return universe, met_pool, gtype.set_axis(eg["gene_id"]).reindex(universe)


@lru_cache(maxsize=3)
def load_inputs(population: str = "primary") -> Inputs:
    cells = population_cells(population)
    lines = cells["CCLE_ID"].tolist()
    universe, met_pool, gtype = _universe()

    sel = pd.read_csv(P1 / "target_selection.csv")
    targets = sorted(sel.loc[sel["target_strict"], "feature"])
    diagnostics = sorted(sel.loc[sel["pipeline_diagnostic"], "feature"])
    classes = sel.set_index("feature")["class"]

    metab = load.load_metabolomics().set_index("CCLE_ID")
    Y = metab.loc[lines, targets + diagnostics].astype(float)

    _, logx = genemap.expression_genes(lines)           # logx: all genes x population lines (ENSGR dropped)

    members = pd.read_csv(P1 / "gene_set_members.csv")
    members = members[members["expressed"] & members["feature"].isin(targets + diagnostics)]
    sets_p = members[members["version"] == "primary"].groupby("feature")["ccle_gene_id"].apply(sorted).to_dict()
    sets_e = members[members["version"] == "enzymatic"].groupby("feature")["ccle_gene_id"].apply(sorted).to_dict()
    for f in targets + diagnostics:
        sets_e.setdefault(f, [])
        assert set(sets_p[f]) <= set(universe), f

    X = logx.loc[universe, lines].to_numpy(dtype=np.float32).T   # lines x genes
    return Inputs(cells=cells, Y=Y, X=X, genes=pd.Index(universe), gene_type=gtype,
                  targets=targets, diagnostics=diagnostics, sets_primary=sets_p, sets_enzymatic=sets_e,
                  met_pool=met_pool, classes=classes)


def baseline_design(cells: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Intercept + one-hot lineage_collapsed + one-hot medium_category (full dummies; solved by lstsq)."""
    D = pd.get_dummies(cells[["lineage_collapsed", "medium_category"]].astype(str), dtype=float)
    D.insert(0, "intercept", 1.0)
    return D.to_numpy(dtype=float), list(D.columns)


def strata(cells: pd.DataFrame) -> np.ndarray:
    """lineage_collapsed x medium_category stratum code per line (for the permuted-target control)."""
    key = cells["lineage_collapsed"].astype(str) + "|" + cells["medium_category"].astype(str)
    return pd.factorize(key)[0]
