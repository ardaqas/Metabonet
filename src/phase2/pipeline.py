"""Per-target evaluation: two-stage model (research_plan §7), repeated grouped CV (§8), metrics (§9).

Within each training fold, in this order and using training rows only:
  1. OLS baseline on intercept + one-hot lineage + one-hot medium           -> M0 predictions, residuals
  2. gene features standardised with training-fold mean and SD
  3. ridge of the residuals on the standardised genes, alpha by LOO (ridge.fit)
  4. test prediction = baseline prediction + ridge prediction
The per-fold baseline and per-fold gene statistics depend only on (target, fold) and (fold) respectively,
so they are computed once and shared by M1, M2, every random set and the permuted-target control.
Out-of-fold R^2 is pooled over the 5 folds of a repeat, then averaged over repeats.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import ridge
from .data import Inputs, baseline_design


@dataclass
class FoldStats:
    """Per (repeat, fold) gene statistics on the training rows (shared across targets)."""
    mu: list[list[np.ndarray]]
    sd: list[list[np.ndarray]]


@dataclass
class TargetFolds:
    """Per (repeat, fold) baseline quantities for one outcome vector."""
    base_te: list[list[np.ndarray]]      # M0 prediction on the test rows
    res_tr: list[list[np.ndarray]]       # training residuals y_tr - M0_tr
    r2_m0: np.ndarray                    # pooled out-of-fold R^2 of M0 per repeat
    rho_m0: np.ndarray


@dataclass
class SetResult:
    r2: np.ndarray            # per repeat, pooled OOF R^2 of M0 + genes
    rho: np.ndarray           # per repeat, Spearman of OOF predictions vs y
    alpha_idx: np.ndarray     # (repeats x folds) chosen alpha index
    k: int = 0

    def delta(self, tf: TargetFolds) -> np.ndarray:
        return self.r2 - tf.r2_m0


class Evaluator:
    def __init__(self, inp: Inputs, folds: list[list[tuple[np.ndarray, np.ndarray]]], alphas=ridge.ALPHAS):
        self.inp, self.folds, self.alphas = inp, folds, alphas
        self.D, self.design_cols = baseline_design(inp.cells)
        self.n = len(inp.cells)
        X = inp.X
        mu, sd = [], []
        for reps in folds:
            mu_r, sd_r = [], []
            for tr, _ in reps:
                Xtr = X[tr].astype(np.float64)
                m = Xtr.mean(axis=0)
                s = Xtr.std(axis=0, ddof=0)
                s[s == 0] = 1.0                      # constant gene in a training fold -> zero feature
                mu_r.append(m)
                sd_r.append(s)
            mu.append(mu_r)
            sd.append(sd_r)
        self.stats = FoldStats(mu, sd)

    # ---- baseline --------------------------------------------------------------------------------
    def target_folds(self, y: np.ndarray) -> TargetFolds:
        base_te, res_tr, r2, rho = [], [], [], []
        sst = ((y - y.mean()) ** 2).sum()
        for reps in self.folds:
            b_r, r_r = [], []
            pred = np.empty(self.n)
            for tr, te in reps:
                beta, *_ = np.linalg.lstsq(self.D[tr], y[tr], rcond=None)
                p_te = self.D[te] @ beta
                b_r.append(p_te)
                r_r.append(y[tr] - self.D[tr] @ beta)
                pred[te] = p_te
            base_te.append(b_r)
            res_tr.append(r_r)
            r2.append(1.0 - ((y - pred) ** 2).sum() / sst)
            rho.append(spearmanr(pred, y).statistic)
        return TargetFolds(base_te, res_tr, np.array(r2), np.array(rho))

    # ---- gene sets -------------------------------------------------------------------------------
    def _z(self, r: int, f: int, rows: np.ndarray, gidx: np.ndarray) -> np.ndarray:
        X = self.inp.X[np.ix_(rows, gidx)].astype(np.float64)
        return (X - self.stats.mu[r][f][gidx]) / self.stats.sd[r][f][gidx]

    def evaluate(self, y: np.ndarray, tf: TargetFolds, gidx: np.ndarray) -> SetResult:
        sst = ((y - y.mean()) ** 2).sum()
        r2, rho = [], []
        aidx = np.zeros((len(self.folds), len(self.folds[0])), dtype=int)
        for r, reps in enumerate(self.folds):
            pred = np.empty(self.n)
            for f, (tr, te) in enumerate(reps):
                Ztr = self._z(r, f, tr, gidx)
                Zte = self._z(r, f, te, gidx)
                coef, b0, a, _ = ridge.fit(Ztr, tf.res_tr[r][f], self.alphas)
                aidx[r, f] = a
                pred[te] = tf.base_te[r][f] + ridge.predict(Zte, coef, b0)
            r2.append(1.0 - ((y - pred) ** 2).sum() / sst)
            rho.append(spearmanr(pred, y).statistic)
        return SetResult(np.array(r2), np.array(rho), aidx, k=len(gidx))

    def null_draws(self, y: np.ndarray, tf: TargetFolds, pool: np.ndarray, k: int, B: int,
                   rng: np.random.Generator, draws: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """B random size-k subsets of `pool` (gene positions), or the given `draws`.
        Returns (delta R^2 per draw x repeat, draws)."""
        if draws is None:
            draws = np.stack([rng.choice(pool, size=k, replace=False) for _ in range(B)])
        d = np.empty((len(draws), len(self.folds)))
        for b in range(B):
            d[b] = self.evaluate(y, tf, np.sort(draws[b])).r2 - tf.r2_m0
        return d, draws


def variance_deciles(X: np.ndarray, pool: np.ndarray) -> np.ndarray:
    """Decile (0-9) of expression variance for every gene, computed within `pool` (outcome-blind)."""
    v = X[:, pool].astype(np.float64).var(axis=0)
    ranks = np.argsort(np.argsort(v))
    dec_pool = np.minimum((ranks * 10) // len(pool), 9)
    dec = np.full(X.shape[1], -1, dtype=int)
    dec[pool] = dec_pool
    return dec


def matched_draws(pool: np.ndarray, dec: np.ndarray, target_genes: np.ndarray, B: int,
                  rng: np.random.Generator) -> np.ndarray:
    """B size-k subsets of `pool` whose variance-decile profile equals that of `target_genes` (D47)."""
    by_dec = {d: pool[dec[pool] == d] for d in range(10)}
    need = {d: int((dec[target_genes] == d).sum()) for d in range(10)}
    out = np.empty((B, len(target_genes)), dtype=int)
    for b in range(B):
        parts = [rng.choice(by_dec[d], size=n, replace=False) for d, n in need.items() if n]
        out[b] = np.sort(np.concatenate(parts))
    return out


def permute_within_strata(y: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    yp = y.copy()
    for s in np.unique(strata):
        i = np.flatnonzero(strata == s)
        yp[i] = y[rng.permutation(i)]
    return yp


def target_seed(name: str, salt: str) -> int:
    """Deterministic 32-bit seed from a feature name and a purpose label."""
    import hashlib
    return int(hashlib.sha256(f"{name}|{salt}".encode()).hexdigest()[:8], 16)
