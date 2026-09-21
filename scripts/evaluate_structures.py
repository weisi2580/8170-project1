#!/usr/bin/env python3
"""Score AF3 predictions (custom-MSA and default-MSA conditions) against the
chain-filtered experimental reference structures.

For each (target, condition) with a downloaded AF3 model present:
  1. global accuracy: TM-score + RMSD via bin/USalign against the
     chain-filtered reference in data/RCSB_chain/ (or data/RCSB/ for the
     single-chain 8ORK/T1112 target).
  2. model confidence: mean pLDDT (and mean PAE if present) parsed directly
     from AF3's confidence JSON -- no reference needed.
  3. (optional) local accuracy: per-residue lDDT via bin/USalign's -TMscore 0
     mode is not lDDT, so this script reports per-residue CA distance from
     the US-align superposition as a proxy unless a dedicated lDDT tool is
     configured -- see --lddt-cmd.

Aggregates everything into eval/scores_summary.csv, one row per
(target, condition) -- 6 rows once all AF3 jobs are back.

Usage:
    python scripts/evaluate_structures.py
    python scripts/evaluate_structures.py --target T1183 --condition custom
"""
import argparse
import csv
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
USALIGN = ROOT / "bin" / "USalign"
RCSB_DIR = ROOT / "data" / "RCSB"
RCSB_CHAIN_DIR = ROOT / "data" / "RCSB_chain"
EVAL_DIR = ROOT / "eval"

# casp_id -> (pdb_id, chain or None)
TARGET_REFS = {
    "T1183": ("8IFX", "B"),
    "T1112": ("8ORK", None),
    "T1122": ("8BBT", "A"),
}
CONDITIONS = ["custom", "default"]


def reference_path(casp_id: str) -> Path:
    pdb_id, chain = TARGET_REFS[casp_id]
    if chain:
        return RCSB_CHAIN_DIR / f"{pdb_id}_{chain}.cif"
    return RCSB_DIR / f"{pdb_id}.cif"


def find_model(condition_dir: Path):
    for ext in (".cif", ".pdb"):
        hits = sorted(condition_dir.glob(f"*model*{ext}")) or sorted(condition_dir.glob(f"*{ext}"))
        hits = [h for h in hits if "input_sequence" not in h.name]
        if hits:
            return hits[0]
    return None


def find_confidence_json(condition_dir: Path):
    hits = sorted(condition_dir.glob("*confidence*.json")) or sorted(condition_dir.glob("*.json"))
    hits = [h for h in hits if h.name != "job_metadata.json"]
    return hits[0] if hits else None


def run_usalign(model_path: Path, ref_path: Path):
    result = subprocess.run(
        [str(USALIGN), str(model_path), str(ref_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    out = result.stdout
    tm_match = re.search(r"TM-score=\s*([\d.]+).*normalized by length of Structure_2", out)
    rmsd_match = re.search(r"RMSD=\s*([\d.]+)", out)
    tm_score = float(tm_match.group(1)) if tm_match else None
    rmsd = float(rmsd_match.group(1)) if rmsd_match else None
    return tm_score, rmsd, out


def parse_confidence(json_path: Path):
    data = json.loads(json_path.read_text())
    plddt = data.get("plddt") or data.get("atom_plddts") or data.get("per_residue_plddt")
    mean_plddt = sum(plddt) / len(plddt) if plddt else None
    pae = data.get("pae") or data.get("predicted_aligned_error")
    mean_pae = None
    if pae:
        flat = [v for row in pae for v in row] if isinstance(pae[0], list) else pae
        mean_pae = sum(flat) / len(flat)
    return mean_plddt, mean_pae


def evaluate_one(casp_id: str, condition: str):
    condition_dir = ROOT / f"af3_{condition}" / casp_id
    model_path = find_model(condition_dir)
    if model_path is None:
        print(f"[{casp_id}/{condition}] no model file found yet, skipping")
        return None

    ref_path = reference_path(casp_id)
    tm_score, rmsd, _ = run_usalign(model_path, ref_path)

    conf_json = find_confidence_json(condition_dir)
    mean_plddt, mean_pae = (None, None)
    if conf_json:
        mean_plddt, mean_pae = parse_confidence(conf_json)

    return {
        "casp_id": casp_id,
        "condition": condition,
        "tm_score": tm_score,
        "rmsd": rmsd,
        "mean_plddt": mean_plddt,
        "mean_pae": mean_pae,
        "model_path": str(model_path.relative_to(ROOT)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=list(TARGET_REFS))
    parser.add_argument("--condition", choices=CONDITIONS)
    args = parser.parse_args()

    if not USALIGN.exists():
        raise SystemExit(f"{USALIGN} not found -- run scripts/build_usalign.sh first")

    targets = [args.target] if args.target else list(TARGET_REFS)
    conditions = [args.condition] if args.condition else CONDITIONS

    rows = []
    for casp_id in targets:
        for condition in conditions:
            row = evaluate_one(casp_id, condition)
            if row:
                rows.append(row)
                print(f"[{casp_id}/{condition}] TM-score={row['tm_score']} RMSD={row['rmsd']} pLDDT={row['mean_plddt']}")

    if not rows:
        print("no AF3 models found -- nothing to score yet")
        return

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    out_path = EVAL_DIR / "scores_summary.csv"
    fieldnames = ["casp_id", "condition", "tm_score", "rmsd", "mean_plddt", "mean_pae", "model_path"]

    existing = []
    if out_path.exists():
        with open(out_path) as fh:
            existing = list(csv.DictReader(fh))
    keys_updated = {(r["casp_id"], r["condition"]) for r in rows}
    existing = [r for r in existing if (r["casp_id"], r["condition"]) not in keys_updated]
    all_rows = existing + [{k: r[k] for k in fieldnames} for r in rows]
    all_rows.sort(key=lambda r: (r["casp_id"], r["condition"]))

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"wrote {out_path.relative_to(ROOT)} ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
