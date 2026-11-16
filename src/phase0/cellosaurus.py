"""Cellosaurus look-ups for problematic-line flags and same-donor relationships.

Queries the public Cellosaurus API (https://api.cellosaurus.org) in batches by
accession (RRID / CVCL_xxxx), caches raw text responses under
data/raw/cellosaurus/, and parses them. Re-running uses the cache unless
``refresh=True``.
"""
from __future__ import annotations

import datetime as dt
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests

from .sources import CELLOSAURUS_API, CELLOSAURUS_DIR

FIELDS = "id,ac,acas,sy,oi,hi,cc,ca,di,sx"
BATCH = 50
USER_AGENT = "metabonet-phase0-audit/0.1 (research data audit; python-requests)"


def release_info() -> dict:
    r = requests.get(f"{CELLOSAURUS_API}/release-info?format=json", timeout=60,
                     headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    return r.json()["Cellosaurus"]["header"]["release"]


def fetch(accessions: list[str], refresh: bool = False, pause: float = 1.0) -> list[Path]:
    """Download Cellosaurus text records for ``accessions`` in batches; return cache paths."""
    CELLOSAURUS_DIR.mkdir(parents=True, exist_ok=True)
    accs = sorted({a for a in accessions if isinstance(a, str) and a.startswith("CVCL_")})
    meta_path = CELLOSAURUS_DIR / "query_meta.json"
    paths = []
    for i in range(0, len(accs), BATCH):
        chunk = accs[i : i + BATCH]
        p = CELLOSAURUS_DIR / f"batch_{i // BATCH:03d}.txt"
        paths.append(p)
        if p.exists() and not refresh and p.read_text().count("\n//") >= 1:
            continue
        q = "ac:(" + " OR ".join(chunk) + ")"
        r = requests.get(
            f"{CELLOSAURUS_API}/search/cell-line",
            params={"q": q, "fields": FIELDS, "rows": len(chunk) + 10, "format": "txt"},
            timeout=120, headers={"User-Agent": USER_AGENT},
        )
        r.raise_for_status()
        p.write_text(r.text)
        time.sleep(pause)
    meta = {
        "api": CELLOSAURUS_API,
        "fields": FIELDS,
        "n_accessions": len(accs),
        "n_batches": len(paths),
        "queried_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    if refresh or not meta_path.exists():
        try:
            meta["release"] = release_info()
        except Exception as e:  # pragma: no cover - network
            meta["release"] = f"unavailable: {e}"
        meta_path.write_text(json.dumps(meta, indent=2))
    return paths


def fetch_single(acc: str) -> str:
    r = requests.get(f"{CELLOSAURUS_API}/cell-line/{acc}", params={"format": "txt", "fields": FIELDS},
                     timeout=60, headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    return r.text


def parse(text: str) -> list[dict]:
    """Parse Cellosaurus flat-file text into records."""
    recs, cur = [], None
    for line in text.splitlines():
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("//"):
            if cur:
                recs.append(cur)
            cur = None
            continue
        tag, val = line[:2], line[5:].strip()
        if tag == "ID":
            cur = {"ID": val, "AC": "", "AS": [], "SY": [], "OI": [], "HI": [], "problematic": [],
                   "CA": "", "DI": [], "SX": ""}
        elif cur is None:
            continue
        elif tag == "AC":
            cur["AC"] = val
        elif tag == "AS":
            cur["AS"] += [v.strip() for v in val.split(";") if v.strip()]
        elif tag == "SY":
            cur["SY"] += [v.strip() for v in val.split(";") if v.strip()]
        elif tag == "OI":
            cur["OI"].append(val.split("!")[0].strip())
        elif tag == "HI":
            cur["HI"].append(val.split("!")[0].strip())
        elif tag == "CC" and val.startswith("Problematic cell line:"):
            cur["problematic"].append(val.replace("Problematic cell line:", "").strip())
        elif tag == "CA":
            cur["CA"] = val
        elif tag == "DI":
            cur["DI"].append(val.split(";")[-1].strip())
        elif tag == "SX":
            cur["SX"] = val
    if cur:
        recs.append(cur)
    return recs


_DERIV_RX = re.compile(r"[Ss]hown to be (?:a |an )?(.+?) derivative")


def load_records(paths: list[Path]) -> pd.DataFrame:
    recs = []
    for p in paths:
        recs += parse(p.read_text())
    df = pd.DataFrame(recs).drop_duplicates("AC")
    df["problematic_text"] = df["problematic"].map(lambda v: " | ".join(v))
    df["is_problematic"] = df["problematic"].map(bool)
    df["problematic_type"] = df["problematic"].map(
        lambda v: ";".join(sorted({t.split(".")[0].strip() for t in v})) if v else ""
    )
    df["derivative_of_name"] = df["problematic_text"].map(
        lambda t: (_DERIV_RX.search(t).group(1) if _DERIV_RX.search(t) else "")
    )
    return df


def accession_index(df: pd.DataFrame) -> dict[str, str]:
    """Map primary and secondary accessions -> primary accession."""
    idx = {a: a for a in df["AC"]}
    for ac, sec in zip(df["AC"], df["AS"]):
        for s in sec:
            idx.setdefault(s, ac)
    return idx


def rrids_for_metabolomics() -> list[str]:
    """RRIDs (from current Model.csv) of all CCLE 2019 metabolomics lines, incl. alias-mapped IDs."""
    from . import load

    metab = load.load_metabolomics()
    anno = load.load_annotations().set_index("CCLE_ID")
    model, _ = load.load_model()
    ids = set(metab["DepMap_ID"].fillna(metab["CCLE_ID"].map(anno["depMapID"])).dropna())
    keep = model["ModelID"].isin(ids)
    if "ModelIDAlias" in model.columns:
        keep |= model["ModelIDAlias"].fillna("").map(lambda s: any(x in ids for x in re.split(r"[;,\s]+", s) if x))
    return sorted(model.loc[keep, "RRID"].dropna().unique())


def main(argv=None) -> int:
    import argparse

    from .download import hash_file, load_manifest, save_manifest, upsert
    from .sources import REPO

    ap = argparse.ArgumentParser(description="Fetch Cellosaurus records for metabolomics lines")
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args(argv)
    rrids = rrids_for_metabolomics()
    print(f"{len(rrids)} RRIDs -> Cellosaurus API")
    paths = fetch(rrids, refresh=args.refresh)
    recs = load_records(paths)
    found = set(accession_index(recs))
    missing = [r for r in rrids if r not in found]
    for acc in missing:  # retry individually (e.g. secondary accessions)
        p = CELLOSAURUS_DIR / f"single_{acc}.txt"
        if not p.exists() or args.refresh:
            p.write_text(fetch_single(acc))
            time.sleep(1.0)
    meta = json.loads((CELLOSAURUS_DIR / "query_meta.json").read_text())
    rel = meta.get("release", {})
    rel_s = f"Cellosaurus {rel.get('version', '?')} ({rel.get('updated', '?')})" if isinstance(rel, dict) else str(rel)
    manifest = load_manifest()
    for p in sorted(CELLOSAURUS_DIR.glob("*.txt")) + [CELLOSAURUS_DIR / "query_meta.json"]:
        sha, md5, n = hash_file(p)
        manifest = upsert(manifest, {
            "key": "cellosaurus", "path": str(p.relative_to(REPO)), "url": f"{CELLOSAURUS_API}/search/cell-line",
            "release": rel_s, "bytes": n, "sha256": sha, "md5": md5, "http_last_modified": "", "http_etag": "",
            "retrieved_utc": meta.get("queried_utc", ""), "note": "API query cache (fields: " + FIELDS + ")",
        })
    save_manifest(manifest)
    print(f"records parsed: {len(recs)}; missing after batch: {len(missing)}; {rel_s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
