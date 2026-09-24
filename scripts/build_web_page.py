#!/usr/bin/env python3
"""Build the interactive results page, docs/index.html.

Fills web/template.html's __DATA__ placeholder with, per target: backbone
coordinates of the experimental structure and both top-ranked AF3 models
(models superposed onto the reference with US-align's TM-score rotation over
the CASP evaluation unit; pLDDT kept in the B-factor column), per-residue
pLDDT / Cα error / MSA coverage, and the summary metrics.

The page draws the 3D view with 3Dmol.js (WebGL), loaded from cdnjs, so it
needs an internet connection but no server: open docs/index.html in any
browser. docs/ can also be served as-is by GitHub Pages.

Run after scripts/plot_comparison.py.

Usage:
    python scripts/build_web_page.py
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_structures import _parser, ca_table, eval_unit, reference_path, run_usalign  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "web" / "template.html"
OUT = ROOT / "docs" / "index.html"
TARGETS = ["T1183", "T1112", "T1122"]
CONDITIONS = ["custom", "default"]
BACKBONE = ("N", "CA", "C", "O")


def backbone_pdb(cif: Path, u=None, t=None, lo=None, hi=None) -> str:
    chain = next(iter(_parser.get_structure("s", str(cif))[0]))
    lines, n = [], 0
    for res in chain:
        het, num, _ = res.id
        if het.strip() and res.get_resname() != "MSE":
            continue
        if lo is not None and not lo <= num <= hi:
            continue
        name = "MET" if res.get_resname() == "MSE" else res.get_resname()
        for a in BACKBONE:
            if a not in res:
                continue
            xyz = res[a].coord.astype(float)
            if u is not None:
                xyz = u @ xyz + t
            n += 1
            lines.append("ATOM  %5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           %s"
                         % (n, a, name, num, *xyz, res[a].bfactor, a[0]))
    return "\n".join(lines) + "\nEND\n"


def main():
    summary = pd.read_csv(ROOT / "eval" / "target_summary_table.csv").set_index("casp_id")
    samples = pd.read_csv(ROOT / "eval" / "scores_all_samples.csv")
    data = {}
    for t in TARGETS:
        lo, hi = eval_unit(t)
        ref_path = reference_path(t)
        ref = ca_table(ref_path)
        d = {"ref": backbone_pdb(ref_path, lo=lo, hi=hi)}
        for c in CONDITIONS:
            model_path = ROOT / f"af3_{c}" / t / f"{t}_model.cif"
            model = ca_table(model_path)
            resnums = [r for r in sorted(ref) if lo <= r <= hi and r in model]
            _, _, u, tr = run_usalign(model, ref, resnums)
            d[c] = backbone_pdb(model_path, u, tr)
            pr = pd.read_csv(ROOT / "eval" / t / f"per_residue_{c}.csv")
            d[f"{c}_res"] = {
                "resnum": pr.resnum.tolist(),
                "plddt": pr.plddt.tolist(),
                "err": [None if pd.isna(x) else round(float(x), 2) for x in pr.ca_error],
            }
            d[f"{c}_tm5"] = samples[(samples.casp_id == t) & (samples.condition == c)] \
                .sort_values("sample").tm_score.tolist()
        d["aa"] = "".join(pd.read_csv(ROOT / "eval" / t / "per_residue_custom.csv").aa)
        cov = pd.read_csv(ROOT / "eval" / t / "msa_column_coverage.csv")
        d["cov"] = {"custom": cov.custom.tolist(), "default": cov.af3_default_all.tolist()}
        d["row"] = {k: (v.item() if hasattr(v, "item") else v) for k, v in summary.loc[t].items()}
        d["templates"] = {
            c: json.loads((ROOT / f"af3_{c}" / t / "job_metadata.json").read_text())["template_pdb_ids"]
            for c in CONDITIONS
        }
        data[t] = d

    html = TEMPLATE.read_text().replace("__DATA__", json.dumps(data, separators=(",", ":")))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html)
    print(f"wrote {OUT.relative_to(ROOT)} ({len(html) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
