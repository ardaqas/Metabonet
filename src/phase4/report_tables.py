"""Render docs/case_studies.md from docs/case_studies.template.md and results/phase4/case_*.csv.

    .venv/bin/python -m src.phase4.report_tables
"""
from __future__ import annotations

import pandas as pd

from ..phase0.sources import REPO

OUT = REPO / "results" / "phase4"
CASES = ["2-deoxycytidine", "creatine", "kynurenine"]


def _md(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    lines += ["| " + " | ".join(str(v).replace("|", "/") for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def gene_table(name: str) -> str:
    g = pd.read_csv(OUT / f"case_{name}_genes.csv")
    ex = g["example_reactions"].fillna("").map(lambda s: s.split(" ; ")[0].split(" [")[0][:64])
    return _md(pd.DataFrame({
        "Gene": g["gene"], "Role": g["roles"], "Coefficient (mean ± SD)": [f"{m:+.4f} ± {s:.4f}" for m, s in zip(g["coef_mean"], g["coef_sd"])],
        "Same sign": g["sign_consistency"].map(lambda v: f"{v:.0%}"),
        "Expression, log2(TPM+1): mean (SD)": [f"{m:.2f} ({s:.2f})" for m, s in zip(g["expr_mean_log2tpm1"], g["expr_sd_log2tpm1"])],
        "A Human-GEM reaction": ex}))


def main() -> None:
    tpl = (REPO / "docs" / "case_studies.template.md").read_text()
    for c in CASES:
        tpl = tpl.replace("{{GENES_" + c + "}}", gene_table(c))
    assert "{{" not in tpl
    (REPO / "docs" / "case_studies.md").write_text(tpl)
    print("wrote docs/case_studies.md")


if __name__ == "__main__":
    main()
