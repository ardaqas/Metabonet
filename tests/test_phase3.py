"""Unit tests for Phase 3 statistics helpers."""
import numpy as np
import pandas as pd

from src.phase3 import stats


def test_bh_matches_known_values():
    p = pd.Series([0.01, 0.04, 0.03, 0.20, np.nan, 0.50])
    q = stats.bh(p)
    assert np.isnan(q[4])
    # n = 5 valid: sorted p = .01,.03,.04,.20,.50 -> ranked .05,.075,.0667,.25,.5 -> monotone from the top
    assert np.allclose(q[[0, 2, 1, 3, 5]].to_numpy(), [0.05, 0.0667, 0.0667, 0.25, 0.5], atol=1e-3)


def test_emp_p_definition():
    d = np.array([0.0, 0.01, 0.02, 0.03])
    assert stats.emp_p(d, 0.025) == (1 + 1) / 5         # one draw >= observed
    assert stats.emp_p(d, 0.1) == 1 / 5                 # none
    assert stats.emp_p(d, -1.0) == 1.0                  # all


def test_h2_classification_rules():
    assert stats.h2_class(0.01, 3) == "supported"
    assert stats.h2_class(0.30, 0) == "not supported"
    assert stats.h2_class(0.30, 2) == "partially supported"
    assert stats.h2_class(0.30, 5).startswith("partially supported (>=3")
    assert stats.h2_class(0.01, 1).startswith("partially supported (global")
    assert stats.h2_class(0.01, 0).startswith("partially supported (global")


def test_null_block_counts_and_classification():
    t = pd.DataFrame({"feature": list("abcdefgh"), "dr2_m1": [0.10, 0.08, 0.05, 0.04, 0.03, 0.025, 0.02, 0.015],
                      "nmet_median": [0.01] * 8, "nmet_p95": [0.035] * 8,
                      "network_specific_nmet": [True, True, True, False, False, False, False, False],
                      "above_p95_nmet": [True, True, True, True, False, False, False, False]})
    b = stats.null_test_block(t, "nmet")
    assert b["n_network_specific"] == 3 and b["n_above_own_p95"] == 4 and b["network_specific_targets"] == ["a", "b", "c"]
    assert np.isclose(b["wilcoxon_p_vs_median_null"], 1 / 256)      # all 8 differences positive
    assert b["wilcoxon_p_vs_median_null"] < 0.05 and b["classification"] == "supported"
    assert b["binomial_p"] < 0.01
