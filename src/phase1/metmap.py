"""Metabolite mapping: 136 CCLE polar (HILIC) features -> compartment-free Human-GEM metabolites.

Inputs
    data/curation/polar_feature_ids.tsv   hand-curated species + anchor IDs (KEGG, else ChEBI; HMDB only
                                          where neither exists) and LC-MS isomer concerns
Verification (cached API look-ups, src/phase1/idlookup.py)
    KEGG: curated KEGG ID must exist and one of its names must match the curated species name
    ChEBI: IDs from the curation + KEGG DBLINKS; expanded by one hop of conjugate acid/base relations
           (protonation state is irrelevant, research_plan §3.4); HMDB IDs taken from ChEBI xrefs
Matching
    a Human-GEM metabolite (metsNoComp) is an ID candidate if it shares any KEGG / ChEBI / HMDB ID
    name candidates (exact normalised name match) are reported but never sufficient on their own
Status (research_plan §3.4 step 4)
    ambiguous   - feature represents >1 species (summed / isomer feature)
    unmapped    - no Human-GEM metabolite by ID or name
    1:1 confident - single species, exactly one ID candidate, no conflicting candidate, no enantiomer
                  of it in Human-GEM, no LC-MS isomer concern, KEGG/ChEBI verification passed
    anything else -> needs_review = True with a *provisional* status and a written reason
"""
from __future__ import annotations

import re

import pandas as pd

from ..phase0.sources import REPO
from . import humangem, idlookup

CURATION = REPO / "data" / "curation" / "polar_feature_ids.tsv"


def _same_formula(a: str, b: str) -> bool:
    """Formulas equal ignoring hydrogen count and charge (protonation state is irrelevant)."""
    strip = lambda f: re.sub(r"H\d*(?![a-z])", "", re.sub(r"[+-]\d*$", "", f or ""))
    return bool(a) and bool(b) and strip(a) == strip(b)


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _split(s: str, n: int) -> list[str]:
    parts = s.split("|") if s else [""] * n
    return parts + [""] * (n - len(parts))


def load_curation() -> pd.DataFrame:
    return pd.read_csv(CURATION, sep="\t", dtype=str).fillna("")


def species_rows(cur: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in cur.itertuples(index=False):
        sp = r.species.split("|")
        n = len(sp)
        for s, k, c, h in zip(sp, _split(r.kegg, n), _split(r.chebi, n), _split(r.hmdb, n)):
            rows.append({"feature": r.feature, "species": s, "kegg_curated": k, "chebi_curated": c,
                         "hmdb_curated": h, "n_species": n, "lcms_concern": r.lcms_concern == "yes",
                         "lcms_note": r.lcms_note, "kegg_name_ok": r.kegg_name_ok,
                         "hg_name_alias": r.hg_name_alias, "review_flag": r.review_flag})
    return pd.DataFrame(rows)


def resolve_ids(sp: pd.DataFrame) -> pd.DataFrame:
    """Verify and expand IDs per species (network look-ups are cached)."""
    kegg = idlookup.kegg_records(sp["kegg_curated"].tolist())
    chebi_seed = set(sp["chebi_curated"]) | {c for k in kegg.values() for c in k["chebi"]}
    ch = idlookup.chebi_records(sorted(x for x in chebi_seed if x))
    # one hop of conjugates + enantiomers for every seed
    extra = {c for r in ch.values() for c in r["conjugates"] + r["enantiomers"]}
    ch.update(idlookup.chebi_records(sorted(extra - set(ch))))
    extra2 = {c for e in extra if e in ch for c in ch[e]["conjugates"]}  # conjugates of enantiomers
    ch.update(idlookup.chebi_records(sorted(extra2 - set(ch))))
    idlookup.write_meta()

    out = []
    for r in sp.to_dict("records"):
        notes, kegg_ids, chebi_ids, hmdb_ids, enan = [], set(), set(), set(), set()
        r.update(kegg_name="", kegg_name_check="", chebi_name="", formula="")
        k = r["kegg_curated"]
        if k:
            rec = kegg[k]
            if not rec["found"]:
                notes.append(f"KEGG {k} not found")
            else:
                kegg_ids.add(k)
                ok = any(_norm(r["species"]) in _norm(n) or _norm(n) in _norm(r["species"]) for n in rec["names"])
                r["kegg_name"] = rec["names"][0] if rec["names"] else ""
                r["kegg_name_check"] = ("ok" if ok else
                                        (f"ok (manual: {r['kegg_name_ok']})" if r["kegg_name_ok"] else
                                         f"name differs: KEGG '{'; '.join(rec['names'][:3])}'"))
                r["formula"] = rec["formula"]
                chebi_ids |= set(rec["chebi"])
        if r["chebi_curated"]:
            chebi_ids.add(r["chebi_curated"])
        seeds = sorted(chebi_ids)
        for c in seeds:
            rec = ch.get(c)
            if rec and rec["found"]:
                chebi_ids |= set(rec["conjugates"]) | set(rec["secondary"])
                kegg_ids |= set(rec["kegg"]) if not k else set()
                hmdb_ids |= set(rec["hmdb"])
                for e in rec["enantiomers"]:
                    enan.add(e)
                    if e in ch and ch[e]["found"]:
                        enan |= set(ch[e]["conjugates"])
            elif rec is not None:
                notes.append(f"{c} not found in ChEBI")
        for c in list(chebi_ids):
            rec = ch.get(c)
            if rec and rec["found"]:
                hmdb_ids |= set(rec["hmdb"])
        if r["chebi_curated"] and ch.get(r["chebi_curated"], {}).get("found"):
            r["chebi_name"] = ch[r["chebi_curated"]]["name"]
            r["formula"] = r["formula"] or ch[r["chebi_curated"]].get("formula", "")
        hmdb_ids = {humangem.norm_hmdb(h) for h in hmdb_ids}
        hc = humangem.norm_hmdb(r["hmdb_curated"]) if r["hmdb_curated"] else ""
        r["hmdb_curated_verified"] = (hc in hmdb_ids) if hc else None
        if hc:
            hmdb_ids.add(hc)
        r["hmdb_unverified"] = [hc] if (hc and not r["hmdb_curated_verified"]) else []
        r.update(kegg_ids=sorted(kegg_ids), chebi_ids=sorted(chebi_ids - enan), hmdb_ids=sorted(hmdb_ids),
                 enantiomer_chebi=sorted(enan), id_notes="; ".join(notes))
        out.append(r)
    return pd.DataFrame(out)


def match(sp: pd.DataFrame, hg: pd.DataFrame) -> pd.DataFrame:
    idx: dict[tuple[str, str], set[str]] = {}
    for r in hg.itertuples(index=False):
        for ns in ("kegg", "chebi", "hmdb"):
            for v in getattr(r, ns):
                idx.setdefault((ns, v), set()).add(r.metsNoComp)
    by_name: dict[str, set[str]] = {}
    for r in hg.itertuples(index=False):
        for n in r.names_all.split("; "):
            by_name.setdefault(_norm(n), set()).add(r.metsNoComp)
    name_of = dict(zip(hg["metsNoComp"], hg["name"]))
    hg_hmdb = dict(zip(hg["metsNoComp"], hg["hmdb"]))
    hg_formula = dict(zip(hg["metsNoComp"], hg["formula"]))

    res = []
    for r in sp.to_dict("records"):
        via: dict[str, set[str]] = {}
        for ns, col in (("kegg", "kegg_ids"), ("chebi", "chebi_ids"), ("hmdb", "hmdb_ids")):
            for v in r[col]:
                for m in idx.get((ns, v), ()):
                    via.setdefault(m, set()).add(ns)
        enan = set()
        for v in r["enantiomer_chebi"]:
            enan |= idx.get(("chebi", v), set())
        enan -= set(via)
        qnames = {r["species"], r["feature"], r["kegg_name"], r["chebi_name"], r["hg_name_alias"]}
        names = set().union(*[by_name.get(_norm(q), set()) for q in qnames if q])
        unver = set(r["hmdb_unverified"])
        r["hg_unverified_only"] = sorted(
            m for m, ns in via.items() if ns == {"hmdb"} and not (set(hg_hmdb[m]) - unver) - set(r["hmdb_ids"]) - unver
            and set(hg_hmdb[m]) & unver and not (set(hg_hmdb[m]) & (set(r["hmdb_ids"]) - unver)))
        allc = sorted(set(via) | names)
        r["formula_check"] = "; ".join(
            f"{m}:{'ok' if _same_formula(r['formula'], hg_formula[m]) else 'DIFFERS ' + str(hg_formula[m]) + ' vs ' + r['formula']}"
            for m in allc) if r["formula"] else "no reference formula"
        r["hg_id_candidates"] = sorted(via)
        r["hg_id_via"] = "; ".join(f"{m}({'/'.join(sorted(v))})" for m, v in sorted(via.items()))
        r["hg_name_candidates"] = sorted(names)
        r["hg_name_only_candidates"] = sorted(names - set(via))
        r["hg_enantiomer_present"] = sorted(enan)
        r["hg_candidate_names"] = "; ".join(f"{m}={name_of[m]}" for m in sorted(set(via) | names | enan))
        res.append(r)
    return pd.DataFrame(res)


def classify(sm: pd.DataFrame) -> pd.DataFrame:
    """Feature-level status from species-level matches."""
    rows = []
    for feat, g in sm.groupby("feature", sort=False):
        g = g.to_dict("records")
        r0 = g[0]
        cands = sorted({m for s in g for m in s["hg_id_candidates"]})
        base = {
            "feature": feat,
            "n_species": r0["n_species"],
            "species": " | ".join(s["species"] for s in g),
            "kegg": " | ".join(",".join(s["kegg_ids"]) for s in g),
            "chebi": " | ".join(",".join(s["chebi_ids"][:6]) for s in g),
            "hmdb": " | ".join(",".join(s["hmdb_ids"][:4]) for s in g),
            "hmdb_source": " | ".join(
                "ChEBI xref" if s["hmdb_ids"] and s["hmdb_curated_verified"] is not False else
                ("curator (unverified)" if s["hmdb_curated"] else "none") for s in g),
            "kegg_name_check": " | ".join(s["kegg_name_check"] or "n/a" for s in g),
            "id_notes": " | ".join(s["id_notes"] for s in g if s["id_notes"]),
            "hg_metsNoComp": ";".join(cands),
            "hg_match_detail": " | ".join(s["hg_id_via"] or "-" for s in g),
            "hg_candidate_names": " | ".join(s["hg_candidate_names"] or "-" for s in g),
            "hg_name_only_candidates": ";".join(sorted({m for s in g for m in s["hg_name_only_candidates"]})),
            "hg_enantiomer_present": ";".join(sorted({m for s in g for m in s["hg_enantiomer_present"]})),
            "formula_check": " | ".join(s["formula_check"] for s in g),
            "review_flag": r0["review_flag"],
            "lcms_concern": r0["lcms_concern"],
            "lcms_note": r0["lcms_note"],
        }
        status, review, reasons = _decide(g, cands)
        # proposed single Human-GEM target: the unique ID candidate, or the unique name-only candidate
        # for a provisional name-only match
        prop = ""
        if status.startswith("1:1"):
            namec = sorted({m for s in g for m in s["hg_name_only_candidates"]})
            prop = cands[0] if len(cands) == 1 else (namec[0] if not cands and len(namec) == 1 else "")
        base.update(status=status, needs_review=review, reason=reasons, hg_target_proposed=prop)
        rows.append(base)
    return pd.DataFrame(rows)


def _decide(g: list[dict], cands: list[str]) -> tuple[str, bool, str]:
    r = g[0]
    if r["n_species"] > 1:
        per = [len(s["hg_id_candidates"]) for s in g]
        detail = ", ".join(f"{s['species']}: {len(s['hg_id_candidates'])} HG match(es)" for s in g)
        return "ambiguous", False, f"summed/multi-species LC-MS feature ({r['n_species']} species; {detail})"
    s = r
    idc, namec, enan = s["hg_id_candidates"], s["hg_name_only_candidates"], s["hg_enantiomer_present"]
    bad_kegg = "" if (s["kegg_name_check"] == "" or s["kegg_name_check"].startswith("ok")) else s["kegg_name_check"]
    if not idc and not namec:
        return "unmapped", False, "no Human-GEM metabolite shares a KEGG/ChEBI/HMDB ID or name"
    if not idc and namec:
        fok = len(namec) == 1 and "DIFFERS" not in s["formula_check"] and s["formula"]
        if fok and not s["review_flag"] and not s["lcms_concern"]:   # D32
            return ("1:1 confident", False,
                    f"name match {namec[0]} (Human-GEM entry has no KEGG/ChEBI/HMDB ID), formula consistent, "
                    "no competing candidate (D32)")
        if fok:
            return ("1:1 confident (provisional)", True,
                    f"name-only match {namec[0]}; formula consistent; " + (s["review_flag"] or s["lcms_note"]))
        return "unmapped", True, f"name-only match {namec}: no shared ID and formula not confirmed; verify manually"
    issues = []
    unv = s["hg_unverified_only"]
    if unv and not (len(unv) == 1 and unv[0] in s["hg_name_candidates"] and "DIFFERS" not in s["formula_check"]
                    and s["formula"] and len(idc) == 1):   # D32: name + formula rescue
        issues.append(f"match {unv} rests only on a curator-supplied HMDB ID that could not be verified")
    if "DIFFERS" in s["formula_check"]:
        issues.append(f"formula check: {s['formula_check']}")
    if s["review_flag"]:
        issues.append(s["review_flag"])
    if len(idc) > 1:
        issues.append(f"{len(idc)} Human-GEM metabolites share IDs ({', '.join(idc)})")
    if namec:
        issues.append(f"additional name-only candidate(s) {namec}")
    if enan:
        issues.append(f"enantiomer present in Human-GEM ({', '.join(enan)}); LC-MS does not resolve enantiomers")
    if s["lcms_concern"]:
        issues.append(f"LC-MS isomer concern: {s['lcms_note']}")
    if bad_kegg:
        issues.append(f"KEGG {bad_kegg}")
    if s["id_notes"]:
        issues.append(s["id_notes"])
    if not issues:
        return "1:1 confident", False, f"single species; unique Human-GEM match {idc[0]} via {s['hg_id_via']}"
    prov = "ambiguous" if len(set(idc) | set(namec)) > 1 else "1:1 confident (provisional)"
    return prov, True, "; ".join(issues)


DECISIONS = REPO / "data" / "curation" / "review_decisions.tsv"


def apply_decisions(f: pd.DataFrame) -> pd.DataFrame:
    """Apply the project owner's review decisions (D31) to the rule-based classification."""
    dec = pd.read_csv(DECISIONS, sep="\t", dtype=str).fillna("").set_index("feature")
    f = f.copy()
    f["status_rule_based"] = f["status"]
    f["review_decision"] = f["feature"].map(dec["decision"]).fillna("")
    f["review_rationale"] = f["feature"].map(dec["rationale"]).fillna("")
    f["review_ref"] = f["feature"].map(dec["ref"]).fillna("")
    open_items = f.loc[f["needs_review"] & ~f["feature"].isin(dec.index), "feature"].tolist()
    assert not open_items, f"review items without a decision: {open_items}"
    for feat, d in dec.iterrows():
        i = f.index[f["feature"] == feat]
        assert len(i) == 1, feat
        i = i[0]
        if d["decision"] in ("accept", "exclude", "keep_ambiguous"):
            f.at[i, "status"] = d["final_status"]
            f.at[i, "needs_review"] = False
            if d["decision"] == "exclude":
                f.at[i, "hg_target_proposed"] = ""
        elif d["decision"] == "moot":
            f.at[i, "needs_review"] = False          # left provisional; cannot become a target anyway
    return f


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    cur = load_curation()
    sp = resolve_ids(species_rows(cur))
    hg = humangem.metabolite_table()
    sm = match(sp, hg)
    return apply_decisions(classify(sm)), sm
