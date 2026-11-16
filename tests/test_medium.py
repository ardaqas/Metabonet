"""Unit tests for culture-medium normalisation rules (src/phase0/medium.py)."""
import pandas as pd
import pytest

from src.phase0 import medium


@pytest.mark.parametrize("text, base", [
    ("RPMI 1640 + 10% FBS", "RPMI-1640"),
    ("RPMI-1640+ 5% FBS", "RPMI-1640"),
    ("90% rPMI 1640 + 10% h. i. FBS", "RPMI-1640"),
    ("80-90% RPHI 1640+10-20% h.i.FBS", "RPMI-1640"),          # typo in source
    ("RPI-1640+10% FBS", "RPMI-1640"),                          # typo in source
    ("DMEM + 10% FBS", "DMEM"),
    ("DEMEM+10%FBS", "DMEM"),                                   # typo in source
    ("90% Dulbecco's MEM (4.5g/L glucose) + 10% h.i. FBS", "DMEM"),
    ("Dulbecco's modified Eagle's medium with 10% fetal bovine serum", "DMEM"),
    ("EMEM + 10% FBS", "EMEM/MEM"),
    ("MEM + 10% FBS", "EMEM/MEM"),
    ("Eagle's minimal essential medium with 15%  fetal calf serum", "EMEM/MEM"),
    ("AMEM + 20% FBS", "alpha-MEM"),
    ("70% alpha-MEM + 20% h. i. FBS", "alpha-MEM"),
    ("DMEM:F12 = 1:1 + 10% FBS", "DMEM:F12"),
    ("DMEM:HAM's F12  1:1  +5% FBS", "DMEM:F12"),
    ("80% mixture of Dulbecco's MEM + Ham's F 12 at 1:1 + 20% h.i. FBS", "DMEM:F12"),
    ("(DMEM: HamF12=1:1) +10%FBS", "DMEM:F12"),
    ("F12 + 10% FBS", "Ham's F12"),
    ("Han's F12 medium 90% with 10% fetal bovine serum", "Ham's F12"),
    ("F12K + 10% FBS", "Ham's F12K"),
    ("F10 + 10% FBS", "Ham's F10"),
    ("HAMSF10+ 10% FBS", "Ham's F10"),
    ("McCoy's 5A + 10% FBS", "McCoy's 5A"),
    ("McCOY'S5A + 10%FBS", "McCoy's 5A"),
    ("Leibovitz's L-15 + 10% FBS", "Leibovitz L-15"),
    ("L15+10%FBS", "Leibovitz L-15"),
    ("IMDM + 20% FBS", "IMDM"),
    ("80% Iscove's MDM + 20% h.i. FBS", "IMDM"),
    ("William's E Medium + 10% FBS", "Williams' E"),
    ("Waymouth MB 7521 medium+ 10% FBS", "Waymouth's MB 752/1"),
    ("ACL-4", "ACL-4"),
    ("RPMI:F12 = 1:1 + 10% FBS", "mixture"),
    ("DMEM:RPMI = 2:1 + 10% FBS", "mixture"),
    ("MCDB 105:Medium 199 = 1:1 + 15% FBS", "mixture"),
    ("40% Dulbecco's MEM +40% Iscove's MDM + 20% h.i. FBS", "mixture"),
    ("90 % Iscove's MDMD or RPMI 1640 + 10% h.i. FBS", "ambiguous"),
    ("RPMI 1640(or DM-160AU) +10%FBS", "ambiguous"),
    ("", "unknown"),
    (None, "unknown"),
    ("serum-free", "unknown"),
])
def test_classify_base(text, base):
    assert medium.classify(text)["medium_base"] == base


def test_mixture_components_sorted():
    assert medium.classify("RPMI:MEM = 1:1 + 20% FBS")["components"] == "EMEM/MEM+RPMI-1640"
    assert medium.classify("RPMI:EMEM = 1:1 + 20% FBS")["components"] == "EMEM/MEM+RPMI-1640"


@pytest.mark.parametrize("text, pct", [
    ("RPMI 1640 + 10% FBS", 10), ("80% RPMI 1640 +20% h.i. FBS", 20),
    ("RPMI 1640  with L-glutamine (300mg/L), 90%; heat inactivated fetal bovine serum (FBS),10%", 10),
    ("Eagle's minimal essential medium with 15%  fetal calf serum", 15), ("ACL-4", None),
])
def test_serum(text, pct):
    assert medium.serum_percent(text) == pct


def test_assign_prefers_li_and_falls_back():
    idx = pd.Index(["A", "B", "C"])
    li = pd.Series(["RPMI 1640 + 10% FBS", None, None], index=idx)
    an = pd.Series(["DMEM+10%FBS", "EMEM+10%FBS", None], index=idx)
    out = medium.assign(li, an)
    assert out["medium_base"].tolist() == ["RPMI-1640", "EMEM/MEM", "unknown"]
    assert out["medium_source"].tolist() == ["li2019_supp1", "anno2019_growth_medium", "none"]
    assert not bool(out["medium_sources_agree"].iloc[0])     # RPMI (Li) vs DMEM (2019) disagree
    assert pd.isna(out["medium_sources_agree"].iloc[1])


def test_collapse_rule():
    base = pd.Series(["RPMI-1640"] * 20 + ["DMEM"] * 14 + ["mixture"] * 16 + ["unknown"])
    in_set = pd.Series([True] * len(base))
    cat, counts = medium.collapse(base, in_set)
    assert set(cat[:20]) == {"RPMI-1640"}
    assert set(cat[20:34]) == {"other"}          # < 15 lines
    assert set(cat[34:50]) == {"other"}          # mixtures always "other"
    assert cat.iloc[-1] == "unknown"
