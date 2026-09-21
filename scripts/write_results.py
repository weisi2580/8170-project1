#!/usr/bin/env python3
"""Generate RESULTS.md from eval/target_summary_table.csv, stating per target
whether the custom MSA improved accuracy over AF3's default alignment, plus
an overall conclusion across the 3 targets.

Run this last, after scripts/plot_comparison.py.

Usage:
    python scripts/write_results.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT / "eval" / "target_summary_table.csv"
OUT_PATH = ROOT / "RESULTS.md"


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def verdict(tm_custom, tm_default):
    if tm_custom is None or tm_default is None:
        return "not yet evaluated (AF3 model and/or reference score missing)"
    delta = tm_custom - tm_default
    if abs(delta) < 0.01:
        return f"no meaningful difference (ΔTM-score = {delta:+.3f})"
    winner = "custom MSA" if delta > 0 else "AF3 default MSA"
    return f"{winner} was more accurate (ΔTM-score = {delta:+.3f})"


def main():
    if not TABLE_PATH.exists():
        OUT_PATH.write_text(
            "# Results\n\n"
            "No evaluation data yet -- run the pipeline through "
            "scripts/evaluate_structures.py and scripts/plot_comparison.py first.\n"
        )
        print(f"wrote placeholder {OUT_PATH.relative_to(ROOT)} (no eval data yet)")
        return

    with open(TABLE_PATH) as fh:
        rows = list(csv.DictReader(fh))

    lines = ["# Results", ""]
    lines.append(
        "Custom MMseqs2/ColabFold MSA vs. AlphaFold3 default alignment, evaluated "
        "against the experimental PDB reference with US-align (TM-score, RMSD) and "
        "AF3's own confidence (mean pLDDT)."
    )
    lines.append("")

    deltas = []
    for row in rows:
        tm_c = fnum(row.get("tm_score_custom"))
        tm_d = fnum(row.get("tm_score_default"))
        lines.append(f"## {row['casp_id']} (PDB {row.get('pdb_id', '?')}, {row.get('casp_difficulty_class', 'difficulty TBD')})")
        lines.append(
            f"- MSA depth / Neff/L: {row.get('msa_depth', '?')} / {row.get('msa_neff_per_l', '?')}"
        )
        lines.append(
            f"- TM-score: custom={row.get('tm_score_custom', '?')}  default={row.get('tm_score_default', '?')}"
        )
        lines.append(
            f"- RMSD (Å): custom={row.get('rmsd_custom', '?')}  default={row.get('rmsd_default', '?')}"
        )
        lines.append(
            f"- mean pLDDT: custom={row.get('plddt_custom', '?')}  default={row.get('plddt_default', '?')}"
        )
        lines.append(f"- **Verdict:** {verdict(tm_c, tm_d)}")
        lines.append("")
        if tm_c is not None and tm_d is not None:
            deltas.append(tm_c - tm_d)

    lines.append("## Overall conclusion")
    if deltas:
        mean_delta = sum(deltas) / len(deltas)
        n_improved = sum(1 for d in deltas if d > 0.01)
        n_worse = sum(1 for d in deltas if d < -0.01)
        n_flat = len(deltas) - n_improved - n_worse
        lines.append(
            f"Across {len(deltas)}/{len(rows)} fully-evaluated targets, the custom MSA "
            f"improved TM-score on {n_improved}, hurt it on {n_worse}, and made no "
            f"meaningful difference on {n_flat} (mean ΔTM-score = {mean_delta:+.3f})."
        )
    else:
        lines.append(
            "Not all targets have both AF3 conditions evaluated yet -- re-run "
            "scripts/write_results.py once all 6 AF3 jobs (3 targets x 2 conditions) "
            "are downloaded and scored."
        )

    OUT_PATH.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
