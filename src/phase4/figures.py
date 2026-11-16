"""Case-study panels (research_plan §14). Reads results/phase4/case_*.csv.

    .venv/bin/python -m src.phase4.figures
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ..phase0.sources import REPO

OUT = REPO / "results" / "phase4"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SURFACE, INK, INK2, GRID, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#9a9890"
ROLE_COLOR = {"producing": BLUE, "consuming": ORANGE, "transport": AQUA, "reversible / mixed": YELLOW}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlecolor": INK, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "font.size": 9, "figure.dpi": 130, "legend.frameon": False,
})


def role_group(roles: str) -> str:
    r = set(roles.split("/"))
    if r == {"producing"}:
        return "producing"
    if r == {"consuming"}:
        return "consuming"
    if "transport" in r and len(r) == 1:
        return "transport"
    return "reversible / mixed"


def case_figure(name: str, why: str, s: dict) -> None:
    g = pd.read_csv(OUT / f"case_{name}_genes.csv")
    o = pd.read_csv(OUT / f"case_{name}_oof_repeat0.csv")
    g["grp"] = g["roles"].map(role_group)
    g = g.iloc[::-1].reset_index(drop=True)                      # largest |coef| at the top
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 0.36 * len(g) + 2.3), gridspec_kw={"width_ratios": [1.15, 1]})
    y = np.arange(len(g))
    a.barh(y, g["coef_mean"], xerr=g["coef_sd"], color=[ROLE_COLOR[x] for x in g["grp"]], height=0.66,
           edgecolor=SURFACE, lw=0.8, error_kw={"ecolor": INK2, "elinewidth": 0.9, "capsize": 2})
    a.axvline(0, color=INK2, lw=0.8)
    a.set_yticks(y)
    a.set_yticklabels([f"{r.gene}  ({r.roles}; {r.sign_consistency:.0%} same sign)" for r in g.itertuples()], fontsize=7.5)
    a.grid(axis="y", visible=False)
    a.set_xlabel("ridge coefficient on the baseline residual\n(per SD of log2(TPM+1); mean ± SD over the 50 training folds)")
    a.set_title(f"{name}: genes and their direction", loc="left", pad=26)
    a.legend(handles=[Line2D([], [], marker="s", color=c, ls="", ms=7, label=r) for r, c in ROLE_COLOR.items()],
             loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=4, fontsize=7.5, columnspacing=1.2, handletextpad=0.3,
             title="role of the gene's reaction(s) relative to the metabolite", title_fontsize=7.5)
    note = f"largest |coefficient| = {g['coef_mean'].abs().max():.3g}"
    if s["alpha_median_idx"] == 16:
        note += "\nmedian-fold penalty at grid maximum (1e6): M1 = baseline"
    a.text(0.0, 1.0, note, transform=a.transAxes, fontsize=7.2, color=INK2, va="bottom", ha="left")
    # right: observed vs out-of-fold predicted (repeat 0)
    lo, hi = min(o["observed"].min(), o["pred_m1"].min()), max(o["observed"].max(), o["pred_m1"].max())
    b.scatter(o["observed"], o["pred_m0"], s=7, color=MUTED, alpha=0.45, lw=0, label=f"M0 baseline, R² {s['r2_m0']:.3f}", zorder=2)
    b.scatter(o["observed"], o["pred_m1"], s=7, color=BLUE, alpha=0.55, lw=0,
              label=f"M1 + genes, R² {s['r2_m0'] + s['dr2_m1']:.3f}", zorder=3)
    b.plot([lo, hi], [lo, hi], color=INK2, lw=0.8, zorder=1)
    b.set_xlabel("observed (log10 relative abundance)")
    b.set_ylabel("out-of-fold prediction")
    b.set_title(f"ΔR² = {s['dr2_m1']:+.3f}; N-met p = {s['p_nmet']:.3g}\n(R² = 10-repeat mean; points = repeat 0)", loc="left", fontsize=9)
    b.legend(loc="lower right", fontsize=7.5, markerscale=2.2)
    fig.suptitle(f"{why}: {name} (k = {s['k']} genes)", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    OUT.mkdir(parents=True, exist_ok=True)
    stem = f"fig6_case_{name.replace(' ', '_')}"
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    summ = json.loads((OUT / "case_studies_summary.json").read_text())
    for name, s in summ.items():
        case_figure(name, s["why"].capitalize(), s)
    print(f"wrote {len(summ)} case-study figures to {OUT}")


if __name__ == "__main__":
    main()
