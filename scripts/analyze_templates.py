#!/usr/bin/env python3
"""How good were the structural templates AF3 used in each condition?

Both conditions ran with templates ON (useStructureTemplate: true); the AF3
Server searched templates per job and returns the (up to 4) hits it used in
af3_<cond>/<ID>/server_output/templates/, with a query->template residue map
(queryIndices / templateIndices, 0-based; templateIndices index the
template's entity sequence, i.e. label_seq_id - 1).

For every template hit this scores the template itself as if it were a
(partial) model of the target, against the experimental structure:
  - n_aligned / coverage: query residues (inside the CASP evaluation unit)
    the template provides coordinates for
  - seq_identity over the aligned residues
  - tm_ref: TM-score normalized by the full evaluation-unit length, using
    US-align's residue-index superposition of the aligned Cα pairs -- "how
    much of the right answer this template alone hands AF3"
  - rmsd over the aligned residues

and, per (target, condition), the union coverage of all templates plus the
coverage of any extra region passed with --region.

Writes eval/template_analysis.csv (one row per template hit) and
eval/template_summary.csv (one row per target x condition).

Usage:
    python scripts/analyze_templates.py
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import index_to_one, three_to_index

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_structures import ca_table, eval_unit, reference_path, run_usalign  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = ROOT / "eval"
TARGETS = ["T1183", "T1112", "T1122"]
CONDITIONS = ["custom", "default"]


def one_letter(name):
    try:
        return index_to_one(three_to_index("MET" if name == "MSE" else name))
    except (KeyError, ValueError):
        return "X"


def template_ca(cif: Path):
    """label_seq_id -> (one-letter aa, CA xyz) for the template's polymer chain,
    first model, first altloc."""
    m = MMCIF2Dict(str(cif))
    seq = "".join(one_letter(x) for x in m["_entity_poly_seq.mon_id"])
    out = {}
    cols = zip(m["_atom_site.label_atom_id"], m["_atom_site.label_seq_id"],
               m["_atom_site.Cartn_x"], m["_atom_site.Cartn_y"], m["_atom_site.Cartn_z"],
               m.get("_atom_site.pdbx_PDB_model_num", ["1"] * len(m["_atom_site.id"])))
    for atom, sid, x, y, z, model in cols:
        if atom != "CA" or sid in (".", "?") or model != "1":
            continue
        sid = int(sid)
        if sid not in out:
            out[sid] = (seq[sid - 1], np.array([float(x), float(y), float(z)]))
    return seq, out


def tm_over_reference(model_xyz, ref_xyz, u, t, L_ref):
    d0 = 1.24 * (max(L_ref, 22) - 15) ** (1 / 3) - 1.8
    d = np.linalg.norm(model_xyz @ u.T + t - ref_xyz, axis=1)
    return float(np.sum(1 / (1 + (d / d0) ** 2)) / L_ref)


def main():
    hit_rows, summary_rows = [], []
    for casp_id in TARGETS:
        lo, hi = eval_unit(casp_id)
        ref = ca_table(reference_path(casp_id))
        unit = [r for r in sorted(ref) if lo <= r <= hi]
        query = "".join(l.strip() for l in (ROOT / "data" / "CASP15" / f"{casp_id}.fasta")
                        .read_text().splitlines() if not l.startswith(">"))
        for cond in CONDITIONS:
            tdir = ROOT / f"af3_{cond}" / casp_id / "server_output" / "templates"
            maps = sorted(tdir.glob("*_query_to_hit.json"))
            if not maps:
                print(f"[{casp_id}/{cond}] no template data, skipping")
                continue
            covered_union = set()
            best_tm = 0.0
            for rank, hit in enumerate(json.loads(maps[0].read_text())):
                cif = tdir / hit["name"]
                pdb_id = cif.read_text().splitlines()[0].removeprefix("data_")
                seq, tca = template_ca(cif)
                pairs = [(qi + 1, ti + 1) for qi, ti in zip(hit["queryIndices"], hit["templateIndices"])]
                ident = np.mean([query[q - 1] == seq[t - 1] for q, t in pairs])
                usable = [(q, t) for q, t in pairs if q in ref and lo <= q <= hi and t in tca]
                covered_union |= {q for q, _ in usable}
                tm_ref = rmsd = None
                if len(usable) >= 5:
                    qs = [q for q, _ in usable]
                    model = {q: (ref[q][0], tca[t][1], 0.0) for q, t in usable}
                    tm_local, rmsd, u, tr = run_usalign(model, ref, qs)
                    tm_ref = tm_over_reference(np.array([tca[t][1] for _, t in usable]),
                                               np.array([ref[q][1] for q in qs]), u, tr, len(unit))
                    best_tm = max(best_tm, tm_ref)
                qs_sorted = sorted(q for q, _ in usable)
                hit_rows.append({
                    "casp_id": casp_id, "condition": cond, "rank": rank, "template": pdb_id,
                    "n_aligned": len(usable),
                    "coverage": round(len(usable) / len(unit), 3),
                    "query_range": f"{qs_sorted[0]}-{qs_sorted[-1]}" if qs_sorted else "",
                    "seq_identity": round(float(ident), 3),
                    "tm_ref": None if tm_ref is None else round(tm_ref, 3),
                    "rmsd_aligned": rmsd,
                })
            summary_rows.append({
                "casp_id": casp_id, "condition": cond,
                "templates": " ".join(r["template"] for r in hit_rows
                                      if r["casp_id"] == casp_id and r["condition"] == cond),
                "union_coverage": round(len(covered_union) / len(unit), 3),
                "best_template_tm_ref": round(best_tm, 3),
            })

    for name, rows in [("template_analysis.csv", hit_rows), ("template_summary.csv", summary_rows)]:
        with open(EVAL_DIR / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote eval/{name}")
    for r in hit_rows:
        print(f"{r['casp_id']} {r['condition']:7s} #{r['rank']} {r['template']}  aligned {r['n_aligned']:3d} "
              f"({r['coverage']:.0%}, res {r['query_range']:>8s})  id {r['seq_identity']:.2f}  "
              f"TM_ref {r['tm_ref']}  RMSD {r['rmsd_aligned']}")
    for r in summary_rows:
        print(r)


if __name__ == "__main__":
    main()
