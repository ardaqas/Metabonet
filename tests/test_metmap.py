"""Unit tests for the metabolite-mapping status rules (src/phase1/metmap._decide); no network access."""
from src.phase1 import metmap


def _species(**kw):
    base = {"n_species": 1, "species": "x", "feature": "x", "hg_id_candidates": [], "hg_name_only_candidates": [],
            "hg_name_candidates": [], "hg_enantiomer_present": [], "lcms_concern": False, "lcms_note": "",
            "kegg_name_check": "ok", "id_notes": "", "hg_unverified_only": [], "formula_check": "MAM1:ok",
            "formula": "C6H12O6", "review_flag": "", "hg_id_via": "MAM1(kegg)"}
    base.update(kw)
    return base


def test_unique_id_match_is_confident():
    s, review, _ = metmap._decide([_species(hg_id_candidates=["MAM1"])], ["MAM1"])
    assert s == "1:1 confident" and not review


def test_d32_name_formula_no_competitor_is_confident():
    sp = _species(hg_id_candidates=[], hg_name_only_candidates=["MAM1"], hg_name_candidates=["MAM1"], hg_id_via="")
    s, review, reason = metmap._decide([sp], [])
    assert s == "1:1 confident" and not review and "D32" in reason


def test_name_match_with_formula_mismatch_is_unmapped_for_review():
    sp = _species(hg_id_candidates=[], hg_name_only_candidates=["MAM1"], hg_name_candidates=["MAM1"],
                  formula_check="MAM1:DIFFERS C5H10O5 vs C6H12O6", hg_id_via="")
    s, review, _ = metmap._decide([sp], [])
    assert s == "unmapped" and review


def test_name_match_with_lcms_concern_stays_provisional():
    sp = _species(hg_id_candidates=[], hg_name_only_candidates=["MAM1"], hg_name_candidates=["MAM1"],
                  lcms_concern=True, lcms_note="isobaric", hg_id_via="")
    s, review, _ = metmap._decide([sp], [])
    assert s == "1:1 confident (provisional)" and review


def test_unverified_hmdb_rescued_by_name_and_formula():
    sp = _species(hg_id_candidates=["MAM1"], hg_name_candidates=["MAM1"], hg_unverified_only=["MAM1"],
                  hg_id_via="MAM1(hmdb)")
    s, review, _ = metmap._decide([sp], ["MAM1"])
    assert s == "1:1 confident" and not review


def test_unverified_hmdb_without_name_match_needs_review():
    sp = _species(hg_id_candidates=["MAM1"], hg_unverified_only=["MAM1"], hg_id_via="MAM1(hmdb)")
    s, review, _ = metmap._decide([sp], ["MAM1"])
    assert review


def test_enantiomer_present_needs_review():
    sp = _species(hg_id_candidates=["MAM1"], hg_enantiomer_present=["MAM2"])
    s, review, reason = metmap._decide([sp], ["MAM1"])
    assert s == "1:1 confident (provisional)" and review and "enantiomer" in reason


def test_multi_species_is_ambiguous():
    g = [_species(n_species=2, species="a", hg_id_candidates=["MAM1"]),
         _species(n_species=2, species="b", hg_id_candidates=["MAM2"])]
    s, review, _ = metmap._decide(g, ["MAM1", "MAM2"])
    assert s == "ambiguous" and not review


def test_two_id_candidates_is_ambiguous_for_review():
    sp = _species(hg_id_candidates=["MAM1", "MAM2"], hg_id_via="MAM1(kegg); MAM2(kegg)")
    s, review, _ = metmap._decide([sp], ["MAM1", "MAM2"])
    assert s == "ambiguous" and review
