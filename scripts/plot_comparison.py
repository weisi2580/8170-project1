#!/usr/bin/env python3
"""Step 7 comparison figures + final target summary table.

Reads eval/scores_summary.csv (from scripts/evaluate_structures.py) and
data/metadata.csv, and writes:
  - figures/metric_comparison.png   (TM-score & RMSD, custom vs. default, per target)
  - figures/<ID>/plddt.png          (per-residue pLDDT, both conditions overlaid),
                                     only if AF3 confidence JSONs are present
  - eval/target_summary_table.csv   (all metrics + CASP difficulty class per target)

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
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = ROOT / "eval"
FIG_DIR = ROOT / "figures"

sns.set_theme(style="whitegrid")


def load_csv(path: Path):
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def metric_comparison(scores):
    if not scores:
        print("eval/scores_summary.csv is empty -- run scripts/evaluate_structures.py first")
        return
    targets = sorted({r["casp_id"] for r in scores})
    by_key = {(r["casp_id"], r["condition"]): r for r in scores}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.arange(len(targets))
    width = 0.35

    for ax, metric, ylabel in [
        (axes[0], "tm_score", "TM-score"),
        (axes[1], "rmsd", "RMSD (Å)"),
    ]:
        custom_vals = [float(by_key[(t, "custom")][metric]) if (t, "custom") in by_key and by_key[(t, "custom")][metric] else np.nan for t in targets]
        default_vals = [float(by_key[(t, "default")][metric]) if (t, "default") in by_key and by_key[(t, "default")][metric] else np.nan for t in targets]
        ax.bar(x - width / 2, custom_vals, width, label="Custom MSA")
        ax.bar(x + width / 2, default_vals, width, label="AF3 default")
        ax.set_xticks(x)
        ax.set_xticklabels(targets)
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.legend()

    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "metric_comparison.png", dpi=150)
    plt.close(fig)
    print(f"wrote {(FIG_DIR / 'metric_comparison.png').relative_to(ROOT)}")


def plddt_plots():
    for casp_id_dir in sorted((ROOT / "af3_custom").glob("*")):
        casp_id = casp_id_dir.name
        curves = {}
        for condition in ("custom", "default"):
            cond_dir = ROOT / f"af3_{condition}" / casp_id
            jsons = [f for f in cond_dir.glob("*confidence*.json")] if cond_dir.exists() else []
            if not jsons:
                continue
            data = json.loads(jsons[0].read_text())
            plddt = data.get("plddt") or data.get("atom_plddts") or data.get("per_residue_plddt")
            if plddt:
                curves[condition] = plddt

        if not curves:
            continue

        plt.figure(figsize=(9, 3.5))
        for condition, plddt in curves.items():
            plt.plot(range(1, len(plddt) + 1), plddt, label=f"{condition} MSA")
        plt.xlabel("Residue")
        plt.ylabel("pLDDT")
        plt.ylim(0, 100)
        plt.title(f"{casp_id}: per-residue pLDDT")
        plt.legend()
        plt.tight_layout()
        out_dir = FIG_DIR / casp_id
        out_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_dir / "plddt.png", dpi=150)
        plt.close()
        print(f"wrote {(out_dir / 'plddt.png').relative_to(ROOT)}")


def target_summary_table(scores, metadata, msa_qc):
    meta_by_id = {r["casp_id"]: r for r in metadata}
    qc_by_id = {r["casp_id"]: r for r in msa_qc}
    scores_by_key = {(r["casp_id"], r["condition"]): r for r in scores}

    rows = []
    all_ids = sorted(set(meta_by_id) | {r["casp_id"] for r in scores})
    for casp_id in all_ids:
        meta = meta_by_id.get(casp_id, {})
        qc = qc_by_id.get(casp_id, {})
        custom = scores_by_key.get((casp_id, "custom"), {})
        default = scores_by_key.get((casp_id, "default"), {})
        rows.append(
            {
                "casp_id": casp_id,
                "pdb_id": meta.get("pdb_id", ""),
                "length": meta.get("length", ""),
                "casp_difficulty_class": meta.get("casp_difficulty_class", ""),
                "msa_depth": qc.get("depth", ""),
                "msa_neff_per_l": qc.get("neff_per_l", ""),
                "tm_score_custom": custom.get("tm_score", ""),
                "tm_score_default": default.get("tm_score", ""),
                "rmsd_custom": custom.get("rmsd", ""),
                "rmsd_default": default.get("rmsd", ""),
                "plddt_custom": custom.get("mean_plddt", ""),
                "plddt_default": default.get("mean_plddt", ""),
            }
        )

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    out_path = EVAL_DIR / "target_summary_table.csv"
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_path.relative_to(ROOT)}")


def main():
    scores = load_csv(EVAL_DIR / "scores_summary.csv")
    metadata = load_csv(ROOT / "data" / "metadata.csv")
    msa_qc = load_csv(EVAL_DIR / "msa_qc_summary.csv")

    metric_comparison(scores)
    plddt_plots()
    target_summary_table(scores, metadata, msa_qc)


if __name__ == "__main__":
    main()
