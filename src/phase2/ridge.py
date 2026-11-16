"""Closed-form ridge regression with exact leave-one-out selection of the penalty (D41).

Equivalent to sklearn.linear_model.RidgeCV(alphas, fit_intercept=True, cv=None) — an unpenalised
intercept plus L2-penalised coefficients, with alpha chosen by efficient LOO on the training data —
but implemented with one SVD (p <= n) or one eigendecomposition of the Gram matrix (p > n) per fit so
that ~100,000 fits per target are affordable. Verified against RidgeCV in tests/test_phase2.py.

Inputs must be the *training fold only*: Z (n x p) with columns already standardised on training
statistics, and y (n,). Nothing from the test fold enters.
"""
from __future__ import annotations

import numpy as np

ALPHAS = np.logspace(-2, 6, 17)


def fit(Z: np.ndarray, y: np.ndarray, alphas: np.ndarray = ALPHAS) -> tuple[np.ndarray, float, int, np.ndarray]:
    """Return (coef (p,), intercept, index of chosen alpha, LOO mean squared error per alpha)."""
    n, p = Z.shape
    ybar = float(y.mean())
    yc = y - ybar
    if p <= n:
        U, s, Vt = np.linalg.svd(Z, full_matrices=False)
        s2 = s ** 2
        Uty = U.T @ yc
        # per alpha: fitted = U diag(s2/(s2+a)) U'y ; hat diag = 1/n + sum_j U_ij^2 s2_j/(s2_j+a)
        w = s2[:, None] / (s2[:, None] + alphas[None, :])                   # r x A
        fitted = U @ (w * Uty[:, None])                                       # n x A
        h = 1.0 / n + (U ** 2) @ w                                            # n x A
        loo = ((yc[:, None] - fitted) / (1.0 - h)) ** 2
        loo_mse = loo.mean(axis=0)
        a = int(np.argmin(loo_mse))
        coef = Vt.T @ ((s / (s2 + alphas[a])) * Uty)
    else:
        K = Z @ Z.T
        lam, Q = np.linalg.eigh(K)
        lam = np.clip(lam, 0.0, None)
        Qty = Q.T @ yc
        w = lam[:, None] / (lam[:, None] + alphas[None, :])
        fitted = Q @ (w * Qty[:, None])
        h = 1.0 / n + (Q ** 2) @ w
        loo = ((yc[:, None] - fitted) / (1.0 - h)) ** 2
        loo_mse = loo.mean(axis=0)
        a = int(np.argmin(loo_mse))
        dual = Q @ (Qty / (lam + alphas[a]))                                  # (K + aI)^-1 yc
        coef = Z.T @ dual
    return coef, ybar, a, loo_mse


def predict(Zte: np.ndarray, coef: np.ndarray, intercept: float) -> np.ndarray:
    return intercept + Zte @ coef
