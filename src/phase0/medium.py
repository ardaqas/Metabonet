"""Culture-medium normalisation (research_plan.md §4 step 5).

Free-text medium descriptions are parsed into base-medium *components* with an
ordered list of regular expressions. Each rule consumes the text it matches so
that, e.g., "Dulbecco's MEM" is not also counted as "MEM".

Classification of a string:
    * contains " or "                       -> "ambiguous"   (not assigned)
    * no component found                    -> "unknown"     (not assigned)
    * {DMEM, Ham's F12}                     -> "DMEM:F12"
    * exactly one component                 -> that component
    * two or more components                -> "mixture"     (components listed)

Base categories (``medium_base``) are deliberately fine-grained; collapsing of
rare categories into "other" (``medium_category``) is done separately with
:func:`collapse` using the < 15-lines rule on the primary analysis set.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

UNASSIGNED = ("unknown", "ambiguous")
MIN_CATEGORY_N = 15


@dataclass(frozen=True)
class Rule:
    component: str
    pattern: str

    @property
    def rx(self) -> re.Pattern:
        return re.compile(self.pattern, re.IGNORECASE)


# Order matters: more specific patterns first; each match is removed from the text.
RULES: list[Rule] = [
    Rule("Ham's F12K", r"f\s*-?\s*12\s*k"),
    Rule("alpha-MEM", r"alpha\s*-?\s*mem|\bamem\b"),
    Rule("DMEM", r"\bd[e]?mem\b|dulbecco'?s?\s*(modified\s*eagle'?s?\s*medium|mem)?|\bdm-160au\b"),
    Rule("IMDM", r"\bimdm\b|iscove'?s?\s*mdmd?"),
    Rule("Ham's F12", r"(ham'?s?|han'?s|hams)?\s*\bf\s*[:\-]?\s*12\b|hamf12|hams\s*f12"),
    Rule("Ham's F10", r"(ham'?s?|hams)?\s*\bf\s*-?\s*10\b|hamf10|hamsf10"),
    Rule("RPMI-1640", r"\brp[mh]?i\b|\brp[mh]i\s*-?\s*1640|\brpi\s*-?\s*1640|\brpm\b|rpmi1640"),
    Rule("McCoy's 5A", r"mccoy'?s?\s*5a"),
    Rule("Leibovitz L-15", r"leibovitz'?s?\s*l\s*-?\s*15|\bl\s*-?\s*15\b"),
    Rule("Williams' E", r"william'?s?\s*e"),
    Rule("Waymouth's MB 752/1", r"waymouth'?s?"),
    Rule("MCDB 105", r"mcdb\s*105"),
    Rule("Medium 199", r"(medium|media)\s*199"),
    Rule("ACL-4", r"\bacl\s*-?\s*4\b"),
    Rule("EMEM/MEM", r"\bemem\b|\bmem\b|eagle'?s?(\s*minimal\s*essential(\s*medium)?)?"),
]

_SERUM_RX = [
    # "... fetal bovine serum (FBS),10%"
    re.compile(r"\((?:fbs|fcs)\)\s*,\s*(\d+(?:\.\d+)?)\s*%", re.I),
    # "10% FBS", "+10 %FBS", "10% h.i. FBS", "10-20% h. i. FBS", "with 10% fetal bovine serum", "+10 FBS"
    re.compile(
        r"(\d+(?:\.\d+)?)(?:\s*-\s*\d+(?:\.\d+)?)?\s*%?\s*(?:h\s*[.,]?\s*i\s*[.,]?\s*)?"
        r"(?:heat[- ]?inactivated\s*)?(?:fbs|fcs|fetal|calf)",
        re.I,
    ),
]


def components(text: str) -> list[str]:
    """Return base-medium components found in ``text`` (in rule order)."""
    if not isinstance(text, str) or not text.strip():
        return []
    s = text
    found: list[str] = []
    for rule in RULES:
        if rule.rx.search(s):
            found.append(rule.component)
            s = rule.rx.sub(" ", s)
    return found


def serum_percent(text: str) -> float | None:
    """First serum percentage mentioned (informational only)."""
    if not isinstance(text, str):
        return None
    for rx in _SERUM_RX:
        m = rx.search(text)
        if m:
            return float(m.group(1))
    return None


def classify(text: str) -> dict:
    """Normalise one free-text medium string."""
    if not isinstance(text, str) or not text.strip():
        return {"medium_base": "unknown", "components": "", "rule": "missing", "serum_pct": None}
    comps = components(text)
    serum = serum_percent(text)
    if re.search(r"\bor\b", text, re.I) and len(comps) > 1:
        return {"medium_base": "ambiguous", "components": "+".join(comps),
                "rule": "alternatives given ('or')", "serum_pct": serum}
    if not comps:
        return {"medium_base": "unknown", "components": "", "rule": "no base medium recognised",
                "serum_pct": serum}
    if set(comps) == {"DMEM", "Ham's F12"}:
        return {"medium_base": "DMEM:F12", "components": "DMEM+Ham's F12",
                "rule": "DMEM + Ham's F12 -> DMEM:F12", "serum_pct": serum}
    if len(comps) == 1:
        return {"medium_base": comps[0], "components": comps[0], "rule": f"single component: {comps[0]}",
                "serum_pct": serum}
    return {"medium_base": "mixture", "components": "+".join(sorted(comps)),
            "rule": "two or more base media", "serum_pct": serum}


def normalise_series(s: pd.Series) -> pd.DataFrame:
    """Vectorised :func:`classify` over a Series (index preserved)."""
    out = pd.DataFrame([classify(v) for v in s], index=s.index)
    return out


def assign(li_text: pd.Series, anno_text: pd.Series) -> pd.DataFrame:
    """Combine sources: Li et al. Supp. Table 1 preferred, 2019 Growth.Medium as fallback.

    Both Series must share an index (CCLE_ID).
    """
    li = normalise_series(li_text)
    an = normalise_series(anno_text)
    use_li = ~li["medium_base"].isin(UNASSIGNED)
    use_an = ~use_li & ~an["medium_base"].isin(UNASSIGNED)
    out = pd.DataFrame(index=li_text.index)
    out["medium_raw_li2019"] = li_text
    out["medium_raw_anno2019"] = anno_text
    out["medium_base_li2019"] = li["medium_base"]
    out["medium_base_anno2019"] = an["medium_base"]
    out["medium_components"] = li["components"].where(use_li, an["components"].where(use_an, ""))
    out["medium_base"] = li["medium_base"].where(use_li, an["medium_base"].where(use_an, "unknown"))
    out["medium_source"] = "li2019_supp1"
    out.loc[~use_li, "medium_source"] = "anno2019_growth_medium"
    out.loc[~use_li & ~use_an, "medium_source"] = "none"
    out["serum_pct"] = li["serum_pct"].where(use_li, an["serum_pct"])
    out["medium_sources_agree"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
    both = ~li["medium_base"].isin(UNASSIGNED) & ~an["medium_base"].isin(UNASSIGNED)
    agree = (li["medium_base"] == an["medium_base"]) & (
        (li["medium_base"] != "mixture") | (li["components"] == an["components"])
    )
    out.loc[both, "medium_sources_agree"] = agree[both]
    return out


def collapse(base: pd.Series, in_set: pd.Series, min_n: int = MIN_CATEGORY_N) -> tuple[pd.Series, pd.Series]:
    """Collapse base categories with < ``min_n`` lines *within the set* into "other".

    Returns (collapsed categories for all rows, counts of base categories in the set).
    Unassigned rows stay "unknown".
    """
    counts = base[in_set].value_counts()
    keep = set(counts[counts >= min_n].index) - set(UNASSIGNED) - {"mixture"}
    cat = base.where(base.isin(keep), "other")
    cat = cat.where(~base.isin(UNASSIGNED), "unknown")
    return cat, counts


def normalisation_table(li_text: pd.Series, anno_text: pd.Series) -> pd.DataFrame:
    """One row per unique raw string x source with its normalisation."""
    rows = []
    for source, s in (("li2019_supp1", li_text), ("anno2019_growth_medium", anno_text)):
        vc = s.fillna("<missing>").value_counts()
        for raw, n in vc.items():
            c = classify(None if raw == "<missing>" else raw)
            rows.append({"source": source, "raw_string": raw, "n_lines": int(n), **c})
    df = pd.DataFrame(rows)
    return df.sort_values(["source", "medium_base", "n_lines"], ascending=[True, True, False]).reset_index(drop=True)
