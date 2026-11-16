"""Download Phase 0 raw data and record provenance + checksums.

Usage (from repo root):
    .venv/bin/python -m src.phase0.download            # fetch missing files, verify all
    .venv/bin/python -m src.phase0.download --force    # re-download everything
    .venv/bin/python -m src.phase0.download --verify   # only re-hash and compare with MANIFEST

Raw files are untrusted data: they are only ever read as data (pandas/openpyxl),
never executed or imported.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
from pathlib import Path

import pandas as pd
import requests

from .sources import HUMANGEM_BLOBS, MANIFEST, MODEL_RELEASE_FILE, REPO, SOURCES, Source

USER_AGENT = "metabonet-phase0-audit/0.1 (research data audit; python-requests)"
MANIFEST_COLS = [
    "key", "path", "url", "release", "bytes", "sha256", "md5",
    "http_last_modified", "http_etag", "retrieved_utc", "note",
]


def hash_file(path: Path) -> tuple[str, str, int]:
    sha, md5, n = hashlib.sha256(), hashlib.md5(), 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)
            md5.update(chunk)
            n += len(chunk)
    return sha.hexdigest(), md5.hexdigest(), n


def git_blob_sha1(path: Path) -> str:
    """SHA-1 as computed by git for a blob (lets us check files against the pinned commit)."""
    data = path.read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def load_manifest() -> pd.DataFrame:
    if MANIFEST.exists():
        return pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
    return pd.DataFrame(columns=MANIFEST_COLS)


def upsert(manifest: pd.DataFrame, row: dict) -> pd.DataFrame:
    manifest = manifest[manifest["path"] != row["path"]]
    return pd.concat([manifest, pd.DataFrame([row])], ignore_index=True)[MANIFEST_COLS]


def save_manifest(manifest: pd.DataFrame) -> None:
    manifest.sort_values("path").to_csv(MANIFEST, sep="\t", index=False)


def fetch(src: Source) -> dict:
    """Stream ``src.url`` to ``src.path`` atomically; return HTTP provenance."""
    src.path.parent.mkdir(parents=True, exist_ok=True)
    tmp = src.path.with_suffix(src.path.suffix + ".part")
    print(f"  downloading {src.url}")
    with requests.get(src.url, stream=True, timeout=120, headers={"User-Agent": USER_AGENT}) as r:
        r.raise_for_status()
        ctype = r.headers.get("Content-Type", "")
        if "text/html" in ctype:
            raise RuntimeError(f"{src.url} returned HTML (possible block page), not data")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
        # Content-Length is the *encoded* size when the server compresses (e.g. GitHub raw); skip then
        expected = None if r.headers.get("Content-Encoding") else r.headers.get("Content-Length")
        info = {
            "http_last_modified": r.headers.get("Last-Modified", ""),
            "http_etag": r.headers.get("ETag", "").strip('"'),
        }
    if expected and tmp.stat().st_size != int(expected):
        raise RuntimeError(f"size mismatch for {src.filename}: {tmp.stat().st_size} != {expected}")
    tmp.replace(src.path)
    info["retrieved_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return info


def model_release() -> str:
    """Release label for the manual DepMap files (README header preferred; see load.release_evidence)."""
    from .load import model_release as _mr

    return _mr()


def release_problems() -> list[str]:
    """Consistency checks between README.txt header, RELEASE.txt and the Model.csv schema."""
    from .load import release_evidence

    ev, out = release_evidence(), []
    if not ev["readme_release"]:
        out.append("data/raw/depmap/README.txt missing or has no '# DepMap Public <release>' header; "
                   "download it from the same DepMap release page as Model.csv.")
    if not ev["release_txt"]:
        out.append(f"MISSING {MODEL_RELEASE_FILE.relative_to(REPO)} (one line naming the DepMap release).")
    if ev["readme_release"] and ev["release_txt"] and ev["readme_release"] not in ev["release_txt"]:
        out.append(f"RELEASE.txt ('{ev['release_txt']}') disagrees with README.txt header ('{ev['readme_release']}').")
    if ev["schema_match"] is False:
        out.append(f"Model.csv columns ({ev['n_file_columns']}) do not match the README's documented Model.csv "
                   f"schema ({ev['n_documented_columns']}); the files may come from different releases.")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="re-download files that already exist")
    ap.add_argument("--verify", action="store_true", help="do not download; only verify checksums")
    args = ap.parse_args(argv)

    manifest = load_manifest()
    problems: list[str] = []

    for src in SOURCES.values():
        rel = str(src.path.relative_to(REPO))
        prev = manifest[manifest["path"] == rel]
        prev_row = prev.iloc[0].to_dict() if len(prev) else {}
        print(f"[{src.key}] {rel}")

        if src.url is None:  # manual DepMap files (Model.csv, README.txt)
            if not src.path.exists():
                problems.append(
                    f"MISSING manual file {rel}. Download {src.filename} from https://depmap.org/portal/data_page/ "
                    f"(current DepMap Public release) in a browser and save it as {rel}; also write the release "
                    f"name (e.g. 'DepMap Public 26Q3, downloaded YYYY-MM-DD') to "
                    f"{MODEL_RELEASE_FILE.relative_to(REPO)}."
                )
                continue
            release = model_release()
            sha, md5, n = hash_file(src.path)
            if prev_row and prev_row.get("sha256") and prev_row["sha256"] != sha:
                msg = f"CHECKSUM CHANGED for {rel}: manifest {prev_row['sha256'][:12]} vs file {sha[:12]}"
                if args.verify:
                    problems.append(msg)
                else:
                    print("  NOTE: " + msg)
            manifest = upsert(manifest, {
                "key": src.key, "path": rel, "url": "manual: https://depmap.org/portal/data_page/",
                "release": release, "bytes": n, "sha256": sha, "md5": md5,
                "http_last_modified": "", "http_etag": "",
                "retrieved_utc": dt.datetime.fromtimestamp(src.path.stat().st_mtime, dt.timezone.utc)
                .strftime("%Y-%m-%dT%H:%M:%SZ") + " (file mtime)",
                "note": src.note,
            })
            print(f"  ok (manual) sha256={sha[:12]}... release={release}")
            continue

        info = {k: prev_row.get(k, "") for k in ("http_last_modified", "http_etag", "retrieved_utc")}
        if args.verify:
            if not src.path.exists():
                problems.append(f"MISSING {rel}")
                continue
        elif args.force or not src.path.exists():
            info = fetch(src)
        sha, md5, n = hash_file(src.path)
        if src.group == "humangem":
            blob = git_blob_sha1(src.path)
            if blob != HUMANGEM_BLOBS[src.filename]:
                problems.append(f"git blob sha1 mismatch for {rel}: {blob} != {HUMANGEM_BLOBS[src.filename]}")
            else:
                print(f"  git blob sha1 verified ({blob[:12]})")
        if prev_row and prev_row.get("sha256") and prev_row["sha256"] != sha:
            msg = f"CHECKSUM CHANGED for {rel}: manifest {prev_row['sha256'][:12]} vs file {sha[:12]}"
            if args.verify:
                problems.append(msg)
            else:
                print("  WARNING: " + msg)
        manifest = upsert(manifest, {
            "key": src.key, "path": rel, "url": src.url, "release": src.release,
            "bytes": n, "sha256": sha, "md5": md5, **info, "note": src.note,
        })
        print(f"  ok sha256={sha[:12]}... bytes={n:,}")

    problems += release_problems()
    if not args.verify:
        save_manifest(manifest)
        print(f"wrote {MANIFEST.relative_to(REPO)}")
    for p in problems:
        print("ERROR: " + p, file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
