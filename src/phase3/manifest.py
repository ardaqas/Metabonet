"""Content hashes of the per-target result arrays (the .npz files are not committed; see .gitignore).

    .venv/bin/python -m src.phase3.manifest            # (re)write results/phase3/RESULTS_MANIFEST.tsv
    .venv/bin/python -m src.phase3.manifest --verify   # check the arrays on disk against the manifest

The hash covers the numeric arrays only: keys starting with "_" (run timings) are excluded, floats are rounded to
6 decimals (the last digits depend on BLAS thread count, ~1e-16; everything reported is far coarser), integers
(random draws, penalty indices) are exact. Bit-identical regeneration additionally needs the same library
versions and OMP_NUM_THREADS=2 (as in `make phase3`).
"""
from __future__ import annotations

import argparse
import hashlib
import sys

import numpy as np
import pandas as pd

from ..phase0.sources import REPO

OUT = REPO / "results" / "phase3"
ANALYSES = ("main", "sens1", "sens2", "sens3", "sens4")
DECIMALS = 6


def array_hash(path) -> str:
    z = np.load(path)
    return arrays_hash({k: z[k] for k in z.files})


def arrays_hash(arrays: dict) -> str:
    h = hashlib.sha256()
    for k in sorted(k for k in arrays if not k.startswith("_")):
        a = np.asarray(arrays[k])
        if a.dtype.kind == "f":
            a = np.round(a.astype(np.float64), DECIMALS) + 0.0          # +0.0 turns -0.0 into 0.0
        h.update(k.encode()); h.update(str(a.dtype).encode()); h.update(str(a.shape).encode())
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    path = OUT / "RESULTS_MANIFEST.tsv"
    if a.verify:
        man = pd.read_csv(path, sep="\t")
        bad = [r.file for r in man.itertuples() if not (REPO / r.file).exists() or array_hash(REPO / r.file) != r.array_sha256]
        print(f"{len(man) - len(bad)} / {len(man)} arrays match the manifest")
        for b in bad[:10]:
            print("MISMATCH or missing:", b, file=sys.stderr)
        return 1 if bad else 0
    rows = [{"file": str(p.relative_to(REPO)), "array_sha256": array_hash(p), "bytes": p.stat().st_size}
            for d in ANALYSES for p in sorted((OUT / d).glob("*.npz"))]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    print(f"hashed {len(rows)} arrays -> results/phase3/RESULTS_MANIFEST.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
