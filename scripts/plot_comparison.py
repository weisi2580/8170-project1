#!/usr/bin/env python3
"""Step 7 comparison figures + final target summary table.

Reads eval/scores_summary.csv, eval/scores_all_samples.csv,
eval/<ID>/per_residue_<cond>.csv (scripts/evaluate_structures.py),
eval/msa_comparison.csv + eval/<ID>/msa_column_coverage.csv
(scripts/compare_msas.py) and data/metadata.csv, and writes:

  figures/metric_comparison.png      TM-score / Cα-lDDT / RMSD per target and
                                     condition: bar = top-ranked sample,
                                     dots = all 5 samples
  figures/msa_depth_vs_accuracy.png  Neff/L of the MSA each condition used vs.
                                     TM-score (the working hypothesis, in one plot)
  figures/plddt_vs_lddt.png          per-residue confidence vs. accuracy
  figures/<ID>/per_residue.png       MSA column coverage, pLDDT and Cα error
                                     along the sequence, both conditions
  figures/<ID>/pae.png               PAE matrices, both conditions
  eval/target_summary_table.csv      one row per target, everything side by side

Usage:
    python scripts/plot_comparison.py
"""
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = ROOT / "eval"
FIG_DIR = ROOT / "figures"

# ordered by MSA depth (deep -> shallow), the axis of the hypothesis
TARGETS = ["T1183", "T1112", "T1122"]
CONDITIONS = ["custom", "default"]
LABELS = {"custom": "Custom MSA", "default": "AF3 default MSA"}
COLORS = {"custom": "#2a78d6", "default": "#eb6834"}  # categorical slots 1, 2
INK, INK_2, INK_3 = "#0b0b0b", "#52514e", "#8a8984"
GRID, SURFACE = "#e4e3df", "#fcfcfb"
PAE_CMAP = LinearSegmentedColormap.from_list(
    "blue_seq", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.linewidth": 1, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 1, "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False, "text.color": INK, "axes.labelcolor": INK_2,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.titlecolor": INK, "font.size": 9.5,
    "legend.frameon": False, "lines.linewidth": 2, "lines.solid_capstyle": "round",
})


def load_csv(path: Path):
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def pick(df, **eq):
    mask = np.ones(len(df), dtype=bool)
    for k, v in eq.items():
        mask &= (df[k] == v).to_numpy()
    return df[mask]


def save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")


def legend_handles():
    return [plt.Line2D([], [], color=COLORS[c], lw=6, label=LABELS[c]) for c in CONDITIONS]


def metric_comparison(summary, samples):
    metrics = [("tm_score", "TM-score", "higher is better", (0, 1)),
               ("lddt", "Cα-lDDT", "higher is better", (0, 1)),
               ("rmsd", "Cα RMSD (Å)", "lower is better", None)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    x = np.arange(len(TARGETS))
    w = 0.34
    for ax, (col, name, note, ylim) in zip(axes, metrics):
        for j, cond in enumerate(CONDITIONS):
            xs = x + (j - 0.5) * w
            top = [pick(summary, casp_id=t, condition=cond)[col].item() for t in TARGETS]
            ax.bar(xs, top, width=w - 0.04, color=COLORS[cond], alpha=0.9, zorder=2)
            for xi, t, v in zip(xs, TARGETS, top):
                vals = pick(samples, casp_id=t, condition=cond)[col].to_numpy()
                ax.scatter(np.full(len(vals), xi), vals, s=16, color=SURFACE,
                           edgecolor=INK, linewidth=0.8, zorder=3)
                ax.annotate(f"{v:.2f}" if col != "rmsd" else f"{v:.1f}", (xi, max(v, vals.max())),
                            ha="center", va="bottom", fontsize=8, color=INK_2,
                            xytext=(0, 4), textcoords="offset points")
        ax.set_xticks(x, TARGETS)
        ax.set_title(f"{name}  ({note})", loc="left")
        ax.grid(axis="x", visible=False)
        if ylim:
            ax.set_ylim(ylim[0], ylim[1] * 1.08)
            ax.set_yticks(np.linspace(*ylim, 6))
    fig.legend(handles=legend_handles() + [plt.Line2D([], [], marker="o", ls="", mfc=SURFACE,
               mec=INK, label="individual AF3 samples (5 per run)")],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.06))
    fig.text(0.5, -0.04, "Bars = top-ranked sample. Targets ordered deep → shallow MSA "
             "(T1183 TBM-easy, T1112 FM/TBM, T1122 FM).", ha="center", color=INK_3, fontsize=8.5)
    fig.tight_layout()
    save(fig, FIG_DIR / "metric_comparison.png")


def msa_depth_vs_accuracy(summary, msa_cmp):
    source_for = {"custom": "custom", "default": "af3_default_all"}
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for t in TARGETS:
        pts = {}
        for cond in CONDITIONS:
            neff = pick(msa_cmp, casp_id=t, msa_source=source_for[cond])["neff_per_l"].item()
            tm = pick(summary, casp_id=t, condition=cond)["tm_score"].item()
            pts[cond] = (neff, tm)
        (x0, y0), (x1, y1) = pts["custom"], pts["default"]
        ax.plot([x0, x1], [y0, y1], color=INK_3, lw=1, zorder=1)
        for cond, (xv, yv) in pts.items():
            ax.scatter(xv, yv, s=70, color=COLORS[cond], edgecolor=SURFACE, linewidth=2, zorder=3)
        ax.annotate(t, (max(x0, x1), max(y0, y1)), xytext=(8, 4), textcoords="offset points",
                    color=INK, fontsize=9.5, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlabel("MSA Neff / L  (80% identity clusters per residue, log scale)")
    ax.set_ylabel("TM-score (top-ranked model)")
    ax.set_ylim(0.3, 1.02)
    ax.set_title("Deeper MSA ↔ better model, but the custom MSA never beat AF3's own", loc="left")
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=8, mfc=COLORS[c], mec=COLORS[c],
                                  label=LABELS[c]) for c in CONDITIONS], loc="lower right")
    fig.tight_layout()
    save(fig, FIG_DIR / "msa_depth_vs_accuracy.png")


def per_residue_plots(meta):
    for t in TARGETS:
        pr = {c: load_csv(EVAL_DIR / t / f"per_residue_{c}.csv") for c in CONDITIONS}
        cov = load_csv(EVAL_DIR / t / "msa_column_coverage.csv")
        if any(df.empty for df in pr.values()) or cov.empty:
            continue
        lo, hi = map(int, meta.loc[t, "eval_unit"].split("-"))
        L = len(cov)
        fig, axes = plt.subplots(3, 1, figsize=(11, 6.4), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1.2, 1.2]})
        cov_cols = {"custom": "custom", "default": "af3_default_all"}
        for cond in CONDITIONS:
            axes[0].plot(cov.resnum, cov[cov_cols[cond]].clip(lower=0.8), color=COLORS[cond])
            axes[1].plot(pr[cond].resnum, pr[cond].plddt, color=COLORS[cond])
            err = pr[cond].dropna(subset=["ca_error"])
            axes[2].plot(err.resnum, err.ca_error, color=COLORS[cond])
        if cov[list(cov_cols.values())].to_numpy().max() > 20:
            axes[0].set_yscale("log")
        else:  # a near-empty MSA reads better on a linear integer axis
            axes[0].yaxis.get_major_locator().set_params(integer=True)
            axes[0].set_ylim(0, None)
        axes[0].set_ylabel("MSA sequences\ncovering residue")
        axes[1].set_ylabel("pLDDT")
        axes[1].set_ylim(0, 100)
        for y, lab in [(90, "very high"), (70, "confident"), (50, "low")]:
            axes[1].axhline(y, color=INK_3, lw=0.8, ls=(0, (1, 0)), alpha=0.5, zorder=0)
            axes[1].text(2, y - 2, lab, va="top", fontsize=7.5, color=INK_3)
        axes[2].set_ylabel("Cα error after\nsuperposition (Å)")
        axes[2].set_xlabel("Residue (CASP FASTA numbering)")
        for ax in axes:
            for a, b in [(0.5, lo - 0.5), (hi + 0.5, L + 0.5)]:
                if b > a:
                    ax.axvspan(a, b, color=GRID, alpha=0.6, lw=0, zorder=0)
            ax.set_xlim(0.5, L + 0.5)
        axes[0].set_title(f"{t}: MSA coverage, confidence and error along the chain "
                          f"(grey = outside CASP evaluation unit {lo}-{hi})", loc="left")
        axes[0].legend(handles=legend_handles(), loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.12))
        fig.tight_layout()
        save(fig, FIG_DIR / t / "per_residue.png")


def pae_plots():
    for t in TARGETS:
        paths = {c: ROOT / f"af3_{c}" / t / f"{t}_confidence.json" for c in CONDITIONS}
        if not all(p.exists() for p in paths.values()):
            continue
        fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3))
        for ax, cond in zip(axes, CONDITIONS):
            pae = np.array(json.loads(paths[cond].read_text())["pae"])
            im = ax.imshow(pae, cmap=PAE_CMAP, vmin=0, vmax=31.75, origin="upper",
                           extent=(0.5, len(pae) + 0.5, len(pae) + 0.5, 0.5))
            ax.set_title(f"{LABELS[cond]} — mean {pae.mean():.1f} Å", loc="left")
            ax.set_xlabel("Scored residue")
            ax.set_ylabel("Aligned residue")
            ax.grid(False)
        cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
        cb.set_label("Predicted aligned error (Å) — light = confident")
        cb.outline.set_visible(False)
        fig.suptitle(f"{t}: AF3 predicted aligned error (top-ranked model)", x=0.02, ha="left",
                     fontweight="bold", fontsize=11)
        save(fig, FIG_DIR / t / "pae.png")


def plddt_vs_lddt():
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
    for ax, t in zip(axes, TARGETS):
        for cond in CONDITIONS:
            df = load_csv(EVAL_DIR / t / f"per_residue_{cond}.csv").dropna(subset=["lddt"])
            r = np.corrcoef(df.plddt, df.lddt)[0, 1]
            ax.scatter(df.plddt, df.lddt, s=9, color=COLORS[cond], alpha=0.55, lw=0,
                       label=f"{LABELS[cond]} (r = {r:.2f})")
        ax.plot([0, 100], [0, 1], color=INK_3, lw=1, zorder=0)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1)
        ax.set_title(t, loc="left")
        ax.set_xlabel("pLDDT (AF3 confidence)")
        ax.legend(loc="upper left", fontsize=8, markerscale=1.8, handletextpad=0.2)
    axes[0].set_ylabel("Cα-lDDT vs. experiment (accuracy)")
    fig.suptitle("Is AF3's confidence honest? Per-residue pLDDT vs. measured accuracy "
                 "(diagonal = perfectly calibrated)", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save(fig, FIG_DIR / "plddt_vs_lddt.png")


def target_summary_table(summary, meta, msa_cmp):
    rows = []
    for t in TARGETS:
        s = {c: pick(summary, casp_id=t, condition=c).iloc[0] for c in CONDITIONS}
        neff = lambda src: pick(msa_cmp, casp_id=t, msa_source=src)  # noqa: E731
        rows.append({
            "casp_id": t,
            "pdb_id": meta.loc[t, "pdb_id"],
            "length": meta.loc[t, "length"],
            "casp_difficulty_class": meta.loc[t, "casp_difficulty_class"],
            "eval_unit": meta.loc[t, "eval_unit"],
            "msa_depth_custom": neff("custom")["depth"].item(),
            "msa_depth_default": neff("af3_default_all")["depth"].item(),
            "neff_per_l_custom": neff("custom")["neff_per_l"].item(),
            "neff_per_l_default": neff("af3_default_all")["neff_per_l"].item(),
            **{f"{m}_{c}": s[c][m] for m in ("tm_score", "lddt", "rmsd", "mean_plddt") for c in CONDITIONS},
            "tm_score_mean5_custom": s["custom"]["tm_score_mean5"],
            "tm_score_mean5_default": s["default"]["tm_score_mean5"],
            "delta_tm_custom_minus_default": round(s["custom"]["tm_score"] - s["default"]["tm_score"], 4),
            "delta_lddt_custom_minus_default": round(s["custom"]["lddt"] - s["default"]["lddt"], 4),
        })
    out = EVAL_DIR / "target_summary_table.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(ROOT)}")


def main():
    summary = load_csv(EVAL_DIR / "scores_summary.csv")
    samples = load_csv(EVAL_DIR / "scores_all_samples.csv")
    msa_cmp = load_csv(EVAL_DIR / "msa_comparison.csv")
    meta = load_csv(ROOT / "data" / "metadata.csv").set_index("casp_id")
    if summary.empty or msa_cmp.empty:
        raise SystemExit("run scripts/evaluate_structures.py and scripts/compare_msas.py first")

    metric_comparison(summary, samples)
    msa_depth_vs_accuracy(summary, msa_cmp)
    plddt_vs_lddt()
    per_residue_plots(meta)
    pae_plots()
    target_summary_table(summary, meta, msa_cmp)


if __name__ == "__main__":
    main()
