"""Content hashes of the per-target result arrays (the .npz files are not committed; see .gitignore).

    .venv/bin/python -m src.phase3.manifest
The hash covers the array bytes (sorted keys), not the zip container, so it is stable across runs.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from ..phase0.sources import REPO

OUT = REPO / "results" / "phase3"


def array_hash(path) -> str:
    z = np.load(path)
    h = hashlib.sha256()
    for k in sorted(z.files):
        a = z[k]
        h.update(k.encode()); h.update(str(a.dtype).encode()); h.update(str(a.shape).encode()); h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def main() -> None:
    rows = [{"file": str(p.relative_to(REPO)), "array_sha256": array_hash(p), "bytes": p.stat().st_size}
            for d in ("main", "sens1", "sens2", "sens3", "sens4") for p in sorted((OUT / d).glob("*.npz"))]
    pd.DataFrame(rows).to_csv(OUT / "RESULTS_MANIFEST.tsv", sep="\t", index=False, lineterminator="\n")
    print(f"hashed {len(rows)} arrays -> results/phase3/RESULTS_MANIFEST.tsv")


if __name__ == "__main__":
    main()
