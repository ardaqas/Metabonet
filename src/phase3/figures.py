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


def fig5_other_target_null() -> None:
    """POST HOC (D54): network-specific targets against the other-target null (gene sets of other targets)."""
    d = pd.read_csv(OUT / "other_target_null.csv")
    det = pd.read_csv(OUT / "other_target_null_sets.csv")
    d = d[d["network_specific_nmet"]].sort_values("dr2_m1").reset_index(drop=True)
    n = len(d)
    fig, ax = plt.subplots(figsize=(6.8, 0.2 * n + 1.5))
    y = np.arange(n)
    for i, r in d.iterrows():
        v = det.loc[det["feature"] == r.feature, "dr2"].to_numpy()
        ax.plot([v.min(), v.max()], [i, i], color=MUTED, lw=1.0, zorder=1)
        ax.plot([np.median(v), np.percentile(v, 95)], [i, i], color=ORANGE, lw=3.0, solid_capstyle="butt", zorder=2)
        ax.plot(np.percentile(v, 95), i, marker="|", color=ORANGE, ms=8, mew=1.6, zorder=3)
    judged = ~d["too_few_to_judge"].to_numpy()
    ax.scatter(d["dr2_m1"][judged], y[judged], s=30, color=BLUE, edgecolor=SURFACE, lw=0.8, zorder=4)
    ax.scatter(d["dr2_m1"][~judged], y[~judged], s=34, facecolor=SURFACE, edgecolor=BLUE, lw=1.4, zorder=4)
    ax.axvline(0, color=INK2, lw=0.8, zorder=0)
    ax.set_yticks(y); ax.set_yticklabels([f"{f} (n={k})" for f, k in zip(d["feature"], d["n_eligible"])], fontsize=7)
    ax.set_ylim(-0.8, n - 0.2); ax.grid(axis="y", visible=False)
    ax.set_xlabel("ΔR² (out-of-fold, mean of 10 repeats)")
    ax.set_title("Post hoc: network-specific targets against gene sets of other targets", loc="left", fontsize=10)
    ax.legend(handles=[Line2D([], [], marker="o", color=BLUE, ls="", ms=5, label="target's own first-neighbour set"),
                       Line2D([], [], marker="o", mfc=SURFACE, mec=BLUE, ls="", ms=6, label="own set; < 10 comparison sets (not judged)"),
                       Line2D([], [], color=ORANGE, lw=3, label="other-target sets: median to 95th percentile"),
                       Line2D([], [], color=MUTED, lw=1, label="other-target sets: full range")],
              loc="lower right", fontsize=7.5)
    ax.text(0.0, -0.075, "n = number of eligible other-target sets (disjoint from the target's genes, 0.5–2x its size)",
            transform=ax.transAxes, fontsize=7, color=INK2, va="top")
    _save(fig, "fig5_posthoc_other_target_null")


def fig_readme_scatter() -> None:
    """Compact figure for the README: one dot per target, typical random-set gain (x) against own-gene gain (y)."""
    import matplotlib as mpl
    df = pd.read_csv(OUT / "per_target_main.csv")
    t = df[df["role"] == "target"].copy()
    with mpl.rc_context({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "xtick.labelsize": 12,
                         "ytick.labelsize": 12, "legend.fontsize": 11.5}):
        fig, ax = plt.subplots(figsize=(7.0, 6.6))
        lim = 0.33
        ax.plot([-0.01, lim], [-0.01, lim], color=INK2, lw=1.4, zorder=1)
        ax.text(0.205, 0.172, "y = x\n(no better than\nrandom genes)", fontsize=11.5, color=INK2, ha="left", va="top", rotation=0)
        spec = t["network_specific_nmet"].to_numpy()
        pc = t["positive_control"].to_numpy()
        ax.scatter(t["nmet_median"][~spec & ~pc], t["dr2_m1"][~spec & ~pc], s=48, color=SURFACE, edgecolor=MUTED, lw=1.5, zorder=3)
        ax.scatter(t["nmet_median"][spec & ~pc], t["dr2_m1"][spec & ~pc], s=58, color=BLUE, edgecolor=SURFACE, lw=0.8, zorder=4)
        ax.scatter(t["nmet_median"][pc], t["dr2_m1"][pc], s=120, marker="D", color=ORANGE, edgecolor=INK, lw=1.1, zorder=5)
        for r in t.nlargest(5, "dr2_m1").itertuples():
            dx, dy, ha = {"2-deoxycytidine": (0.012, 0.0, "left"), "alpha-glycerophosphocholine": (0.012, 0.0, "left"),
                          "2-deoxyadenosine": (0.012, 0.0, "left"), "creatine": (0.012, -0.004, "left"),
                          "taurine": (0.012, -0.004, "left")}.get(r.feature, (0.012, 0.0, "left"))
            name = {"alpha-glycerophosphocholine": "α-glycerophosphocholine"}.get(r.feature, r.feature)
            ax.annotate(name, (r.nmet_median, r.dr2_m1), xytext=(r.nmet_median + dx, r.dr2_m1 + dy), fontsize=12, color=INK,
                        va="center", ha=ha, zorder=6)
        ax.set_xlim(-0.01, lim); ax.set_ylim(-0.01, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("Typical gain from random metabolic genes\n(median ΔR² of 1,000 random sets)")
        ax.set_ylabel("Gain from the metabolite's own network genes (ΔR²)")
        ax.set_title("74 metabolites: own genes against random genes", loc="left", fontweight="bold")
        n_blue, n_open, n_pc = int((spec & ~pc).sum()), int((~spec & ~pc).sum()), int(pc.sum())
        n_pc_win = int((spec & pc).sum())
        handles = [Line2D([], [], marker="o", color=BLUE, ls="", ms=8, markeredgecolor=SURFACE, label=f"clear win over random genes ({n_blue + n_pc_win}, incl. {n_pc_win} controls)"),
                   Line2D([], [], marker="o", mfc=SURFACE, mec=MUTED, mew=1.5, ls="", ms=8, label=f"not a clear win ({n_open})"),
                   Line2D([], [], marker="D", color=ORANGE, ls="", ms=9, markeredgecolor=INK,
                          label=f"positive controls ({n_pc}; {n_pc_win} clear wins)")]
        ax.legend(handles=handles, loc="lower right", borderaxespad=0.4, handletextpad=0.3, fontsize=10.2)
        FIG.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG / "fig_readme_scatter.png", bbox_inches="tight")
        plt.close(fig)


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
    n_extra = 0
    if (OUT / "other_target_null.csv").exists():          # post hoc, D54
        fig5_other_target_null()
        n_extra = 1
    fig_readme_scatter()
    print(f"wrote {5 + n_extra} figures to {FIG}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "readme":      # only the README scatter (leaves the other figure files untouched)
        fig_readme_scatter()
        print(f"wrote {FIG / 'fig_readme_scatter.png'}")
    else:
        main()
