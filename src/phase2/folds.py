"""Cross-validation folds (research_plan §8; D42).

5 folds x 10 repeats, stratified by lineage_collapsed, grouped by donor_group_id (same-donor lines never
straddle a fold boundary). Folds are generated once, saved to results/phase2/folds.csv, and reused
unchanged by every model, every random set and every sensitivity analysis.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

N_SPLITS, N_REPEATS, SEED0 = 5, 10, 1000


def make_folds(cells: pd.DataFrame, n_splits: int = N_SPLITS, n_repeats: int = N_REPEATS,
               seed0: int = SEED0) -> pd.DataFrame:
    y = cells["lineage_collapsed"].astype(str).to_numpy()
    groups = cells["donor_group_id"].astype(str).to_numpy()
    rows = []
    for r in range(n_repeats):
        sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed0 + r)
        for f, (_, te) in enumerate(sgkf.split(np.zeros(len(cells)), y, groups)):
            rows += [{"repeat": r, "fold": f, "CCLE_ID": cells["CCLE_ID"].iat[i]} for i in te]
    df = pd.DataFrame(rows)
    validate_folds(df, cells, n_splits, n_repeats)
    return df


def validate_folds(folds: pd.DataFrame, cells: pd.DataFrame, n_splits: int = N_SPLITS,
                   n_repeats: int = N_REPEATS) -> None:
    ids = set(cells["CCLE_ID"])
    dg = cells.set_index("CCLE_ID")["donor_group_id"]
    for r, g in folds.groupby("repeat"):
        assert set(g["CCLE_ID"]) == ids and len(g) == len(ids), f"repeat {r}: lines not partitioned"
        assert g["fold"].nunique() == n_splits
        nf = g.assign(dg=g["CCLE_ID"].map(dg)).groupby("dg")["fold"].nunique()
        assert (nf == 1).all(), f"repeat {r}: donor group split across folds"
    assert folds["repeat"].nunique() == n_repeats


def fold_indices(folds: pd.DataFrame, cells: pd.DataFrame) -> list[list[tuple[np.ndarray, np.ndarray]]]:
    """[repeat][fold] -> (train_idx, test_idx) as positions in `cells`."""
    pos = {c: i for i, c in enumerate(cells["CCLE_ID"])}
    out = []
    for r in sorted(folds["repeat"].unique()):
        g = folds[folds["repeat"] == r]
        reps = []
        for f in sorted(g["fold"].unique()):
            te = np.array(sorted(pos[c] for c in g.loc[g["fold"] == f, "CCLE_ID"]), dtype=int)
            tr = np.setdiff1d(np.arange(len(cells)), te)
            reps.append((tr, te))
        out.append(reps)
    return out
