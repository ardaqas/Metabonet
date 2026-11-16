"""Core figures (research_plan §14) from results/phase3/per_target_main.csv and the null draw means.

    .venv/bin/python -m src.phase3.figures
Writes results/phase3/figures/fig{1..4}_*.{png,svg}.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from ..phase0.sources import REPO

OUT = REPO / "results" / "phase3"
FIG = OUT / "figures"

# palette (dataviz reference instance): categorical slots in fixed order; ink tokens for text/chrome
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SURFACE, INK, INK2, GRID, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#9a9890"
CLASS_COLOR = {"synthesised": BLUE, "medium-supplied": ORANGE, "mixed": AQUA}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlecolor": INK, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "font.size": 9, "figure.dpi": 130, "legend.frameon": False,
})


def _save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def _safe(name):
    return name.replace("/", "_").replace(" ", "_")


def fig1_targets_vs_nulls(t: pd.DataFrame, nulls: dict) -> None:
    """Per-target dR2(M1) with the N-met and N-all null distributions (5th-95th percentile + median)."""
    t = t.sort_values("dr2_m1").reset_index(drop=True)
    n = len(t)
    fig, ax = plt.subplots(figsize=(7.2, 0.19 * n + 1.4))
    y = np.arange(n)
    for i, r in t.iterrows():
        for lab, col, off in (("nmet", AQUA, -0.18), ("nall", YELLOW, 0.18)):
            d = nulls[lab][r.feature]
            lo, hi, md = np.percentile(d, 5), np.percentile(d, 95), np.median(d)
            ax.plot([lo, hi], [i + off, i + off], color=col, lw=2.2, solid_capstyle="round", alpha=0.9, zorder=2)
            ax.plot(md, i + off, marker="|", color=col, ms=6, mew=1.4, zorder=3)
    ax.axvline(0, color=INK2, lw=0.8, zorder=1)
    ax.axvline(0.02, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
    spec = t["network_specific_nmet"].to_numpy()
    ax.scatter(t["dr2_m1"][~spec], y[~spec], s=22, color=BLUE, zorder=4, edgecolor=SURFACE, lw=0.8)
    ax.scatter(t["dr2_m1"][spec], y[spec], s=34, color=BLUE, zorder=5, edgecolor=INK, lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{f}{' ●' if pc else ''}" for f, pc in zip(t["feature"], t["positive_control"])], fontsize=7)
    ax.set_ylim(-0.8, n - 0.2)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("ΔR² over lineage + medium baseline (out-of-fold, mean of 10 repeats)")
    ax.set_title("Network gene set (M1) against size-matched random gene sets, per target", loc="left")
    handles = [Line2D([], [], marker="o", color=BLUE, ls="", ms=5, label="M1: first-neighbour genes"),
               Line2D([], [], marker="o", color=BLUE, ls="", ms=6, markeredgecolor=INK, label="network-specific (q < 0.10, ΔR² ≥ 0.02)"),
               Line2D([], [], color=AQUA, lw=2.2, label="N-met null: 5th–95th pct, | median"),
               Line2D([], [], color=YELLOW, lw=2.2, label="N-all null: 5th–95th pct, | median"),
               Line2D([], [], color=MUTED, ls=(0, (3, 3)), label="minimum meaningful effect (0.02)")]
    ax.legend(handles=handles, loc="lower right", fontsize=7.5)
    ax.text(0.0, -0.035, "● = literature positive control (H5)", transform=ax.transAxes, fontsize=7, color=INK2, va="top")
    _save(fig, "fig1_dr2_m1_vs_nulls")


def fig2_m1_vs_m2(t: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    lim = max(t["dr2_m2"].max(), t["dr2_m1"].max()) * 1.08
    ax.plot([0, lim], [0, lim], color=INK2, lw=0.8, zorder=1)
    ax.plot([0, lim], [0, 0.7 * lim], color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
    for c, col in CLASS_COLOR.items():
        g = t[t["class"] == c]
        ax.scatter(g["dr2_m2"], g["dr2_m1"], s=26, color=col, edgecolor=SURFACE, lw=0.8, label=f"{c} (n={len(g)})", zorder=3)
    top = t.nlargest(6, "dr2_m1")
    for r in top.itertuples():
        ax.annotate(r.feature, (r.dr2_m2, r.dr2_m1), xytext=(4, 3), textcoords="offset points", fontsize=7, color=INK2)
    ax.set_xlim(min(0, t["dr2_m2"].min()) - 0.01, lim)
    ax.set_ylim(min(0, t["dr2_m1"].min()) - 0.01, lim)
    ax.set_xlabel("ΔR² all expressed genes (M2, ceiling)")
    ax.set_ylabel("ΔR² first-neighbour genes (M1)")
    ax.set_title("Network set against the all-gene ceiling", loc="left")
    ax.legend(loc="upper left", fontsize=8, title="target class", title_fontsize=8)
    ax.text(lim, 0.7 * lim, " capture ratio 0.7", fontsize=7, color=MUTED, va="center")
    ax.text(lim, lim, " M1 = M2", fontsize=7, color=INK2, va="center")
    _save(fig, "fig2_dr2_m1_vs_m2")


def fig3_baseline(t: pd.DataFrame) -> None:
    t = t.sort_values("r2_m0").reset_index(drop=True)
    n = len(t)
    fig, ax = plt.subplots(figsize=(6.4, 0.17 * n + 1.2))
    ax.barh(np.arange(n), t["r2_m0"], color=[CLASS_COLOR[c] for c in t["class"]], height=0.72, edgecolor=SURFACE, lw=0.8)
    ax.set_yticks(np.arange(n)); ax.set_yticklabels(t["feature"], fontsize=7)
    ax.set_ylim(-0.8, n - 0.2); ax.grid(axis="y", visible=False)
    ax.set_xlabel("R² of lineage + medium baseline (M0, out-of-fold)")
    ax.set_title("How much lineage and medium explain, per target", loc="left")
    ax.legend(handles=[Patch(color=col, label=c) for c, col in CLASS_COLOR.items()], loc="lower right", fontsize=8, title="target class", title_fontsize=8)
    _save(fig, "fig3_r2_m0_by_class")


def fig4_permuted(t: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.hist(t["dr2_perm"], bins=np.arange(-0.02, 0.031, 0.0025), color=BLUE, edgecolor=SURFACE, lw=0.8)
    for x, lab in ((0.02, "0.02 = failure threshold"), (0.0, "")):
        ax.axvline(x, color=MUTED if x else INK2, lw=0.8, ls=(0, (3, 3)) if x else "-")
    ax.text(0.02, ax.get_ylim()[1] * 0.95, " 0.02 (S1 failure threshold)", fontsize=7.5, color=INK2, va="top")
    ax.set_xlabel("ΔR² of M1 on the within-stratum permuted target (P)")
    ax.set_ylabel("targets")
    ax.set_title(f"Permuted-target control: median {t['dr2_perm'].median():+.4f}, max {t['dr2_perm'].max():+.4f}", loc="left")
    _save(fig, "fig4_permuted_control")


def main() -> None:
    df = pd.read_csv(OUT / "per_target_main.csv")
    t = df[df["role"] == "target"].copy()
    z = np.load(OUT / "null_draw_means_main.npz")
    nulls = {}
    for key in z.files:
        lab, name = key.split("__", 1)
        feat = next(f for f in t["feature"] if _safe(f) == name) if any(_safe(f) == name for f in t["feature"]) else None
        if feat:
            nulls.setdefault(lab, {})[feat] = z[key]
    fig1_targets_vs_nulls(t, nulls)
    fig2_m1_vs_m2(t)
    fig3_baseline(t)
    fig4_permuted(t)
    print(f"wrote 4 figures to {FIG}")


if __name__ == "__main__":
    main()
