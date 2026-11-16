"""Unit tests for the Phase 2 pipeline components (no project data needed)."""
import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import RidgeCV

from src.phase2 import folds as fo
from src.phase2 import pipeline as pl
from src.phase2 import ridge
from src.phase2.data import Inputs, baseline_design, strata

rng = np.random.default_rng(7)


@pytest.mark.parametrize("n,p", [(80, 10), (60, 150)])      # p <= n (SVD path) and p > n (Gram path)
def test_ridge_loo_matches_sklearn_ridgecv(n, p):
    Z = rng.normal(size=(n, p))
    Z = (Z - Z.mean(0)) / Z.std(0)
    y = Z[:, :3] @ np.array([1.0, -0.5, 0.25]) + rng.normal(scale=0.8, size=n) + 2.0
    alphas = np.logspace(-2, 4, 13)
    coef, b0, a, loo = ridge.fit(Z, y, alphas)
    ref = RidgeCV(alphas=alphas, fit_intercept=True, cv=None, store_cv_results=True).fit(Z, y)
    assert np.isclose(alphas[a], ref.alpha_)
    assert np.allclose(coef, ref.coef_, atol=1e-6)
    assert np.isclose(b0, ref.intercept_, atol=1e-6)
    # LOO MSE curve identical to sklearn's stored per-alpha errors
    assert np.allclose(loo, ref.cv_results_.mean(axis=0), rtol=1e-6)


def test_ridge_uses_only_training_rows():
    Z = rng.normal(size=(50, 5)); y = rng.normal(size=50)
    c1 = ridge.fit(Z, y)[0]
    c2 = ridge.fit(Z.copy(), y.copy())[0]
    assert np.array_equal(c1, c2)


def _toy_cells(n=120, seed=0):
    r = np.random.default_rng(seed)
    lineage = r.choice(["A", "B", "C"], size=n, p=[0.5, 0.3, 0.2])
    medium = r.choice(["RPMI", "DMEM"], size=n)
    donor = [f"S_{i}" for i in range(n)]
    for i in range(0, 20, 2):           # 10 donor pairs with the same lineage
        donor[i + 1] = donor[i]
        lineage[i + 1] = lineage[i]
    return pd.DataFrame({"CCLE_ID": [f"L{i}" for i in range(n)], "lineage_collapsed": lineage,
                         "medium_category": medium, "donor_group_id": donor})


def test_folds_partition_group_and_stratify():
    cells = _toy_cells()
    f = fo.make_folds(cells, n_splits=5, n_repeats=3, seed0=1)
    fo.validate_folds(f, cells, 5, 3)
    # stratification: each fold holds every lineage
    lin = cells.set_index("CCLE_ID")["lineage_collapsed"]
    per = f.assign(l=f["CCLE_ID"].map(lin)).groupby(["repeat", "fold"])["l"].nunique()
    assert (per == 3).all()
    # repeats differ, but are reproducible
    f2 = fo.make_folds(cells, n_splits=5, n_repeats=3, seed0=1)
    assert f.equals(f2)
    assert not f[f.repeat == 0]["fold"].equals(f[f.repeat == 1]["fold"].reset_index(drop=True))


def test_fold_indices_cover_all_rows():
    cells = _toy_cells()
    f = fo.make_folds(cells, 5, 2, 3)
    idx = fo.fold_indices(f, cells)
    for reps in idx:
        allte = np.concatenate([te for _, te in reps])
        assert sorted(allte) == list(range(len(cells)))
        for tr, te in reps:
            assert not set(tr) & set(te) and len(tr) + len(te) == len(cells)


def _toy_inputs(n=120, p=30, seed=3):
    cells = _toy_cells(n, seed)
    r = np.random.default_rng(seed)
    X = r.normal(size=(n, p)).astype(np.float32)
    genes = pd.Index([f"G{i}" for i in range(p)])
    lin_eff = cells["lineage_collapsed"].map({"A": 0.0, "B": 0.8, "C": -0.5}).to_numpy()
    y = lin_eff + 0.6 * X[:, 0] - 0.4 * X[:, 1] + r.normal(scale=0.5, size=n)
    Y = pd.DataFrame({"t": y}, index=cells["CCLE_ID"])
    return Inputs(cells=cells, Y=Y, X=X, genes=genes, gene_type=pd.Series("", index=genes), targets=["t"],
                  diagnostics=[], sets_primary={"t": ["G0", "G1", "G2"]}, sets_enzymatic={"t": ["G0"]},
                  met_pool=list(genes[:10]), classes=pd.Series({"t": "synthesised"}))


def test_baseline_residuals_mean_zero_and_m0_r2_positive():
    inp = _toy_inputs()
    idx = fo.fold_indices(fo.make_folds(inp.cells, 5, 2, 11), inp.cells)
    ev = pl.Evaluator(inp, idx)
    tf = ev.target_folds(inp.Y["t"].to_numpy())
    for r in range(2):
        for f in range(5):
            assert abs(tf.res_tr[r][f].mean()) < 1e-9        # OLS with intercept: residuals sum to zero
    assert (tf.r2_m0 > 0.05).all()          # toy lineage effect is modest; must be clearly positive


def test_gene_standardisation_uses_training_rows_only():
    inp = _toy_inputs()
    idx = fo.fold_indices(fo.make_folds(inp.cells, 5, 1, 11), inp.cells)
    ev = pl.Evaluator(inp, idx)
    tr, te = idx[0][0]
    assert np.allclose(ev.stats.mu[0][0], inp.X[tr].astype(float).mean(0))
    # changing test rows must not change training statistics or the training design matrix
    inp2 = _toy_inputs()
    inp2.X[te] += 100.0
    ev2 = pl.Evaluator(inp2, idx)
    assert np.allclose(ev2.stats.mu[0][0], ev.stats.mu[0][0]) and np.allclose(ev2.stats.sd[0][0], ev.stats.sd[0][0])
    assert np.allclose(ev2._z(0, 0, tr, np.array([0, 1])), ev._z(0, 0, tr, np.array([0, 1])))


def test_true_genes_beat_baseline_and_random_sets():
    inp = _toy_inputs()
    idx = fo.fold_indices(fo.make_folds(inp.cells, 5, 2, 11), inp.cells)
    ev = pl.Evaluator(inp, idx)
    y = inp.Y["t"].to_numpy()
    tf = ev.target_folds(y)
    m1 = ev.evaluate(y, tf, inp.idx(["G0", "G1", "G2"]))
    assert (m1.delta(tf) > 0.1).all()
    pool = np.arange(3, 30)
    d, draws = ev.null_draws(y, tf, pool, 3, 20, np.random.default_rng(0))
    assert d.shape == (20, 2) and draws.shape == (20, 3)
    assert (d.mean(axis=1) < m1.delta(tf).mean()).all()


def test_permutation_within_strata_preserves_stratum_values():
    inp = _toy_inputs()
    s = strata(inp.cells)
    y = inp.Y["t"].to_numpy()
    yp = pl.permute_within_strata(y, s, np.random.default_rng(1))
    assert not np.array_equal(y, yp)
    for k in np.unique(s):
        assert np.allclose(np.sort(y[s == k]), np.sort(yp[s == k]))


def test_design_has_intercept_and_all_levels():
    cells = _toy_cells()
    D, cols = baseline_design(cells)
    assert cols[0] == "intercept" and D.shape[1] == 1 + 3 + 2
    assert np.allclose(D[:, 1:4].sum(1), 1) and np.allclose(D[:, 4:].sum(1), 1)


def test_target_seed_is_deterministic():
    assert pl.target_seed("x", "nmet") == pl.target_seed("x", "nmet") != pl.target_seed("x", "nall")
