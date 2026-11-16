"""Cached look-ups of KEGG COMPOUND and ChEBI records used to verify curated metabolite IDs.

* KEGG REST:  https://rest.kegg.jp/get/<C-ids>   (names, formula, ChEBI DBLINKS)
* ChEBI API:  https://www.ebi.ac.uk/chebi/backend/api/public/compound/<CHEBI:n>/
  (name, conjugate acid/base and enantiomer relations, KEGG/HMDB cross-references)

HMDB itself (hmdb.ca) returns a Cloudflare challenge to scripted clients, so HMDB IDs are
taken from ChEBI cross-references where available; curator-supplied HMDB IDs are marked unverified.

Raw responses are cached under data/raw/idlookup/ and recorded in data/raw/MANIFEST.tsv.
"""
from __future__ import annotations

import datetime as dt
import json
import re
import time
from pathlib import Path

import requests

from ..phase0.sources import RAW

CACHE = RAW / "idlookup"
UA = "metabonet-phase1-idcheck/0.1 (research metabolite ID verification; python-requests)"
KEGG = "https://rest.kegg.jp/get/"
CHEBI = "https://www.ebi.ac.uk/chebi/backend/api/public/compound/"


def _get(url: str, tries: int = 4, timeout: int = 90) -> requests.Response:
    for i in range(tries):
        try:
            r = requests.get(url, timeout=timeout, headers={"User-Agent": UA})
            if r.status_code in (200, 404):
                return r
        except requests.RequestException:
            pass
        time.sleep(3 * (i + 1))
    raise RuntimeError(f"failed: {url}")


# --------------------------------------------------------------------------- KEGG
def kegg_records(ids: list[str], pause: float = 0.4) -> dict[str, dict]:
    d = CACHE / "kegg"
    d.mkdir(parents=True, exist_ok=True)
    ids = sorted({i for i in ids if i})
    todo = [i for i in ids if not (d / f"{i}.txt").exists()]
    for k in range(0, len(todo), 10):  # KEGG allows up to 10 entries per request
        chunk = todo[k : k + 10]
        r = _get(KEGG + "+".join(f"cpd:{c}" for c in chunk))
        for block in r.text.split("///"):
            m = re.search(r"^ENTRY\s+(C\d{5})", block, re.M)
            if m:
                (d / f"{m.group(1)}.txt").write_text(block.strip() + "\n///\n")
        for c in chunk:  # remember misses
            if not (d / f"{c}.txt").exists():
                (d / f"{c}.txt").write_text("NOT FOUND\n")
        time.sleep(pause)
    return {i: parse_kegg((d / f"{i}.txt").read_text()) for i in ids}


def parse_kegg(text: str) -> dict:
    if text.startswith("NOT FOUND"):
        return {"found": False, "names": [], "formula": "", "chebi": []}
    names, formula, chebi, field = [], "", [], None
    for line in text.splitlines():
        if line[:12].strip():
            field = line[:12].strip()
        val = line[12:].strip()
        if field == "NAME":
            names.append(val.rstrip(";"))
        elif field == "FORMULA":
            formula = val
        elif field == "DBLINKS" and val.startswith("ChEBI:"):
            chebi += [f"CHEBI:{x}" for x in val.split(":", 1)[1].split()]
    return {"found": True, "names": names, "formula": formula, "chebi": chebi}


# --------------------------------------------------------------------------- ChEBI
def chebi_records(ids: list[str], pause: float = 0.2) -> dict[str, dict]:
    d = CACHE / "chebi"
    d.mkdir(parents=True, exist_ok=True)
    out = {}
    for cid in sorted({i for i in ids if i}):
        p = d / f"{cid.replace(':', '_')}.json"
        if not p.exists():
            r = _get(CHEBI + cid + "/")
            p.write_text(r.text if r.status_code == 200 else json.dumps({"not_found": True}))
            time.sleep(pause)
        out[cid] = parse_chebi(json.loads(p.read_text()), cid)
    return out


def _clean(name: str) -> str:
    return re.sub(r"<[^>]+>", "", name or "")


def parse_chebi(d: dict, cid: str) -> dict:
    if d.get("not_found"):
        return {"found": False, "name": "", "formula": "", "conjugates": [], "enantiomers": [], "kegg": [], "hmdb": [],
                "secondary": []}
    rel = d.get("ontology_relations") or {}
    conj, enan = set(), set()
    me = int(cid.split(":")[1])
    for r in (rel.get("incoming_relations") or []) + (rel.get("outgoing_relations") or []):
        other = r["init_id"] if r["final_id"] == me else r["final_id"]
        if r["relation_type"] in ("is conjugate acid of", "is conjugate base of"):
            conj.add(f"CHEBI:{other}")
        elif r["relation_type"] == "is enantiomer of":
            enan.add(f"CHEBI:{other}")
    kegg, hmdb = set(), set()
    for group in (d.get("database_accessions") or {}).values():
        for x in group:
            src, acc = (x.get("source_name") or ""), x.get("accession_number") or ""
            if src.startswith("KEGG COMPOUND") and re.fullmatch(r"C\d{5}", acc):
                kegg.add(acc)
            if src.upper().startswith("HMDB") and acc.startswith("HMDB"):
                hmdb.add(acc)
    cd = d.get("chemical_data") or {}
    formula = cd.get("formula", "") if isinstance(cd, dict) else ""
    return {"found": True, "name": _clean(d.get("ascii_name") or d.get("name")), "formula": formula or "",
            "conjugates": sorted(conj),
            "enantiomers": sorted(enan), "kegg": sorted(kegg), "hmdb": sorted(hmdb),
            "secondary": [f"CHEBI:{s}" if not str(s).startswith("CHEBI") else s for s in (d.get("secondary_ids") or [])]}


def write_meta() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "query_meta.json").write_text(json.dumps({
        "kegg_api": KEGG, "chebi_api": CHEBI,
        "queried_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "note": "HMDB (hmdb.ca) not queryable by script (Cloudflare challenge, HTTP 403).",
    }, indent=2))


def register_manifest() -> None:
    """Record the cached KEGG/ChEBI responses (one row per file) in data/raw/MANIFEST.tsv."""
    from ..phase0.download import hash_file, load_manifest, save_manifest, upsert
    from ..phase0.sources import REPO

    meta = json.loads((CACHE / "query_meta.json").read_text())
    man = load_manifest()
    for p in sorted(CACHE.rglob("*")):
        if not p.is_file():
            continue
        sha, md5, n = hash_file(p)
        src = "kegg" if p.parent.name == "kegg" else ("chebi" if p.parent.name == "chebi" else "meta")
        man = upsert(man, {
            "key": f"idlookup_{src}", "path": str(p.relative_to(REPO)),
            "url": KEGG if src == "kegg" else (CHEBI if src == "chebi" else ""),
            "release": "live API (KEGG REST / ChEBI backend API)", "bytes": n, "sha256": sha, "md5": md5,
            "http_last_modified": "", "http_etag": "", "retrieved_utc": meta.get("queried_utc", ""),
            "note": "cached ID-verification response",
        })
    save_manifest(man)
