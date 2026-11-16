"""Documentation and repository checks.

    python scripts/check_docs.py links     # every relative link/image in README.md and docs/ resolves; README has no absolute URLs
    python scripts/check_docs.py paths     # no absolute personal paths, no file > 50 MB, large data/results excluded by .gitignore
    python scripts/check_docs.py numbers   # every number quoted in README.md is recomputed from the result files
    python scripts/check_docs.py all

Exit status is non-zero if any check fails. Checks that need the raw data (marked) are skipped if it is absent.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
HTML_SRC = re.compile(r"""(?:src|href)=["']([^"']+)["']""")
ABS_PATH = re.compile(r"(/Users/[A-Za-z0-9_.-]+|/private/tmp|/home/[A-Za-z0-9_.-]+|[A-Za-z]:\\Users\\)")
MAX_BYTES = 50 * 1024 * 1024


def tracked() -> list[Path]:
    out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return [ROOT / f for f in out.split("\0") if f]


# --------------------------------------------------------------------------------------------- links
def check_links() -> list[str]:
    errors = []
    files = [README] + sorted((ROOT / "docs").rglob("*.md"))
    n = 0
    for f in files:
        text = f.read_text()
        targets = LINK.findall(text) + HTML_SRC.findall(text)
        for t in targets:
            if t.startswith(("mailto:", "#")):
                continue
            if t.startswith(("http://", "https://")):
                if f == README:
                    errors.append(f"{f.relative_to(ROOT)}: absolute URL not allowed in README: {t}")
                continue
            n += 1
            path = (f.parent / t.split("#")[0].split("?")[0]).resolve()
            if not path.exists():
                errors.append(f"{f.relative_to(ROOT)}: broken link -> {t}")
    print(f"links: {n} relative links checked in {len(files)} files; {len(errors)} problem(s)")
    return errors


# --------------------------------------------------------------------------------------------- paths / hygiene
def check_paths() -> list[str]:
    errors = []
    files = tracked()
    for f in files:
        if f.stat().st_size > MAX_BYTES:
            errors.append(f"{f.relative_to(ROOT)}: {f.stat().st_size / 1e6:.0f} MB > 50 MB")
        try:
            text = f.read_text()
        except (UnicodeDecodeError, OSError):
            continue
        m = ABS_PATH.search(text)
        if m:
            errors.append(f"{f.relative_to(ROOT)}: absolute personal path '{m.group(0)}'")
    biggest = max(files, key=lambda p: p.stat().st_size)
    must_ignore = ["data/raw/ccle2019/CCLE_metabolomics_20190502.csv", "data/raw/humangem/Human-GEM.xml", "data/interim/x.pkl",
                   "results/phase3/main/creatine.npz", "results/phase3/sens4/creatine.npz", ".venv/bin/python"]
    for p in must_ignore:
        r = subprocess.run(["git", "check-ignore", "-q", p], cwd=ROOT)
        if r.returncode != 0:
            errors.append(f".gitignore does not exclude {p}")
    print(f"paths: {len(files)} tracked files; largest {biggest.relative_to(ROOT)} ({biggest.stat().st_size / 1e6:.1f} MB); "
          f"{len(errors)} problem(s)")
    return errors


# --------------------------------------------------------------------------------------------- numbers
def facts() -> dict:
    import pandas as pd
    f: dict = {}
    g = json.loads((ROOT / "results/phase3/global_tests.json").read_text())
    pt = pd.read_csv(ROOT / "results/phase3/per_target_main.csv")
    tg = pt[pt["role"] == "target"]
    cells = pd.read_csv(ROOT / "data/processed/cell_lines_phase0.csv")
    oth = json.loads((ROOT / "results/phase3/other_target_null_summary.json").read_text())
    f["n_targets"] = len(tg)
    f["n_lines"] = int(cells["in_primary_set"].sum())
    f["n_network_specific"] = g["H2"]["n_network_specific"]
    f["n_judged"] = oth["network_specific_judged"]
    f["n_judged_exceed"] = oth["network_specific_exceed_p95"]
    sens = g["sensitivity"]["robustness"]
    f["sens_robust"] = bool(sens["robust"])
    f["n_sens"] = int(sens["n_analyses_available"])
    f["h1_median"] = g["H1"]["median_dr2_m1"]
    f["h1_n_ge"] = g["H1"]["n_dr2_ge_0.02"]
    f["h1_frac"] = g["H1"]["frac_dr2_ge_0.02"]
    f["h1_p"] = g["H1"]["wilcoxon_p"]
    f["h2_above_p95"] = g["H2"]["n_above_own_p95"]
    f["h2_expected"] = g["H2"]["expected_above_p95"]
    f["frac_network_specific"] = f["n_network_specific"] / f["n_targets"]
    k = tg[tg["feature"] == "kynurenine"].iloc[0]
    f["kyn_dr2"] = float(k["dr2_m1"]); f["kyn_p"] = float(k["p_nmet"])
    f["n_metabolites"] = len(pd.read_csv(ROOT / "data/processed/metabolite_inventory.csv"))
    fo = pd.read_csv(ROOT / "results/phase2/folds.csv")
    f["n_folds"] = int(fo["fold"].nunique()); f["n_repeats"] = int(fo["repeat"].nunique())
    f["fold_hidden_frac"] = 1.0 / f["n_folds"]
    exp = (ROOT / "docs/expectations.md").read_text()
    f["prior_text_present"] = "10–25%" in exp
    out = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q", "tests"], cwd=ROOT, capture_output=True, text=True).stdout
    m = re.search(r"(\d+) tests? collected", out)
    f["n_tests"] = int(m.group(1)) if m else -1
    return f


def raw_facts() -> dict | None:
    """Facts that need the raw expression data (skipped if absent)."""
    if not (ROOT / "data/raw/ccle2019/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz").exists():
        return None
    sys.path.insert(0, str(ROOT))
    import numpy as np
    import pandas as pd
    from sklearn.linear_model import LinearRegression, RidgeCV
    from sklearn.preprocessing import StandardScaler
    from src.phase1 import genemap
    from src.phase2 import ridge
    from src.phase2.data import baseline_design, load_inputs

    out: dict = {}
    cells = pd.read_csv(ROOT / "data/processed/cell_lines_phase0.csv")
    lines = cells.loc[cells["in_primary_set"], "CCLE_ID"].tolist()
    gtf = genemap.read_gtf_genes()
    _, logx = genemap.expression_genes(lines)
    for sym in ("IDO1", "TDO2"):
        gid = gtf.loc[gtf["gene_name"] == sym, "gene_id"].iloc[0]
        out[f"{sym}_frac_ge1"] = float((logx.loc[gid].to_numpy() >= 1).mean())
    # independent re-implementation of M1 (scikit-learn only) for creatine, repeats 0-2
    inp = load_inputs("primary")
    fo = pd.read_csv(ROOT / "results/phase2/folds.csv")
    pos = {c: i for i, c in enumerate(inp.cells.CCLE_ID)}
    D, _ = baseline_design(inp.cells)
    name = "creatine"
    y = inp.Y[name].to_numpy(float)
    gidx = inp.idx(inp.sets_primary[name])
    z = np.load(ROOT / "results/phase3/main/creatine.npz")
    ref = (z["r2_m1"] - z["r2_m0"])[:3]
    mine = []
    for r in range(3):
        p1, p0 = np.empty(len(y)), np.empty(len(y))
        gr = fo[fo["repeat"] == r]
        for k in sorted(gr["fold"].unique()):
            te = np.array(sorted(pos[c] for c in gr.loc[gr["fold"] == k, "CCLE_ID"]))
            tr = np.setdiff1d(np.arange(len(y)), te)
            base = LinearRegression(fit_intercept=False).fit(D[tr], y[tr])
            sc = StandardScaler().fit(inp.X[np.ix_(tr, gidx)].astype(float))
            rg = RidgeCV(alphas=ridge.ALPHAS).fit(sc.transform(inp.X[np.ix_(tr, gidx)].astype(float)), y[tr] - base.predict(D[tr]))
            p0[te] = base.predict(D[te])
            p1[te] = p0[te] + rg.predict(sc.transform(inp.X[np.ix_(te, gidx)].astype(float)))
        sst = ((y - y.mean()) ** 2).sum()
        mine.append((1 - ((y - p1) ** 2).sum() / sst) - (1 - ((y - p0) ** 2).sum() / sst))
    out["indep_max_abs_diff"] = float(np.abs(np.array(mine) - ref).max())
    return out


def check_numbers() -> list[str]:
    text = README.read_text().replace("−", "-")
    f = facts()
    rf = raw_facts()
    rows: list[tuple[str, str, str]] = []     # (claim, snippet that must appear in README, value recomputed from files)

    def claim(label: str, snippet: str, computed: str) -> None:
        rows.append((label, snippet, computed))

    claim("targets tested / network-specific", "37 of 74 pre-specified", f"{f['n_network_specific']} of {f['n_targets']} pre-specified")
    claim("cell lines", "735 CCLE solid-tumour cell lines", f"{f['n_lines']} CCLE solid-tumour cell lines")
    claim("sensitivity analyses (robust)", "four sensitivity analyses",
          ("four" if f["n_sens"] == 4 else str(f["n_sens"])) + " sensitivity analyses" if f["sens_robust"] else "NOT ROBUST")
    claim("other-target check", "(33 of 33 judged)", f"({f['n_judged_exceed']} of {f['n_judged']} judged)")
    claim("pipeline diagram: targets", "74 target metabolites", f"{f['n_targets']} target metabolites")
    claim("pipeline diagram: cell lines", "735 solid-tumour cell lines", f"{f['n_lines']} solid-tumour cell lines")
    claim("H1 median gain", "median gain in held-out R² is 0.029", f"median gain in held-out R² is {f['h1_median']:.3f}")
    claim("H1 targets above 0.02", "51 of 74 metabolites (69%)", f"{f['h1_n_ge']} of {f['n_targets']} metabolites ({f['h1_frac']:.0%})")
    claim("H1 p-value", "p = 7e-14", f"p = {f['h1_p']:.0e}".replace("e-0", "e-"))
    claim("H2 network-specific", "37 of 74 metabolites clearly beat", f"{f['n_network_specific']} of {f['n_targets']} metabolites clearly beat")
    claim("H2 above 95th percentile", "39 exceed the 95th percentile", f"{f['h2_above_p95']} exceed the 95th percentile")
    claim("H2 expected by chance", "about 3.7 would by chance", f"about {f['h2_expected']:.1f} would by chance")
    claim("prior written in advance", "10–25% of metabolites beforehand", "10–25% of metabolites beforehand" if f["prior_text_present"] else "(10–25% not found in expectations.md)")
    claim("observed share", "and found 50%", f"and found {f['frac_network_specific']:.0%}")
    claim("kynurenine ΔR² and p", "ΔR² = -0.001, p = 0.87", f"ΔR² = {f['kyn_dr2']:.3f}, p = {f['kyn_p']:.2f}")
    claim("metabolites measured", "74 of the 225 measured metabolites", f"{f['n_targets']} of the {f['n_metabolites']} measured metabolites")
    claim("cross-validation scheme", "5 folds × 10 repeats", f"{f['n_folds']} folds × {f['n_repeats']} repeats")
    claim("fraction hidden per fold", "hiding a fifth of the cell lines", "hiding a fifth of the cell lines" if abs(f["fold_hidden_frac"] - 0.2) < 1e-9 else f"hiding {f['fold_hidden_frac']:.0%}")
    claim("unit tests (folder map)", "89 unit tests", f"{f['n_tests']} unit tests")
    claim("unit tests (reproduce)", "# 89 unit tests", f"# {f['n_tests']} unit tests")
    if rf is None:
        print("numbers: raw data absent; IDO1/TDO2 percentages and the independent re-implementation check are SKIPPED")
    else:
        claim("IDO1 / TDO2 expressed (log2 TPM+1 >= 1)", "18% and 30% of the cell lines",
              f"{rf['IDO1_frac_ge1']:.0%} and {rf['TDO2_frac_ge1']:.0%} of the cell lines")
        claim("independent re-implementation", "agrees to about 1e-15",
              "agrees to about 1e-15" if rf["indep_max_abs_diff"] < 1e-14 else f"max |diff| = {rf['indep_max_abs_diff']:.1e}")
    errors = []
    width = max(len(r[0]) for r in rows)
    for label, snippet, computed in rows:
        in_text = snippet in text
        match = snippet == computed
        status = "ok" if (in_text and match) else ("MISMATCH" if in_text else "NOT IN README")
        if status != "ok":
            errors.append(f"{label}: README says '{snippet}', files give '{computed}' ({status})")
        print(f"  {status:13s} {label:{width}s} | README: {snippet!r}  files: {computed!r}")
    # the one-line result must be present verbatim
    one_line = ("For 37 of 74 pre-specified polar metabolites in 735 CCLE solid-tumour cell lines, genes in the metabolite's "
                "Human-GEM reaction neighbourhood improved held-out prediction beyond lineage and culture medium, more than "
                "size-matched random metabolic genes (pre-registered; robust in four sensitivity analyses). A post hoc check "
                "showed these sets also beat non-overlapping gene sets from other metabolites (33 of 33 judged).")
    if one_line not in " ".join(text.split()):
        errors.append("the 'Result in one line' text is not present verbatim")
    print(f"numbers: {len(rows)} claims checked; {len(errors)} problem(s)")
    return errors


def main() -> int:
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    errors: list[str] = []
    if what in ("links", "all"):
        errors += check_links()
    if what in ("paths", "all"):
        errors += check_paths()
    if what in ("numbers", "all"):
        errors += check_numbers()
    for e in errors:
        print("ERROR:", e, file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
