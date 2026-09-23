#!/usr/bin/env python3
"""Score AF3 predictions (custom-MSA and default-MSA conditions) against the
chain-filtered experimental reference structures.

Run scripts/import_af3_results.py first. For each (target, condition), every
one of the 5 AF3 samples in af3_<cond>/<ID>/server_output/ is scored, so the
comparison is not hostage to a single seed/sample:

  1. global accuracy: TM-score + RMSD via bin/USalign in -TMscore 1 mode
     (residue-index correspondence, i.e. "how good is this model of this
     sequence", not a free structural alignment). TM-score is normalized by
     the reference length. Model and reference share the CASP FASTA
     numbering, and scoring is restricted to the CASP15 evaluation unit
     (data/metadata.csv: eval_unit).
  2. local accuracy: superposition-free C-alpha lDDT (inclusion radius 15 A,
     thresholds 0.5/1/2/4 A) computed here in numpy, per residue and as the
     mean over residues. Cα-only, so it is not identical to the all-atom
     lDDT CASP reports, but it ranks models the same way for this purpose.
  3. model confidence: per-residue pLDDT (C-alpha B-factor of the AF3 model),
     mean PAE from the full_data JSON, pTM + ranking_score from the summary.

Writes:
  eval/scores_all_samples.csv        one row per (target, condition, sample)
  eval/scores_summary.csv            one row per (target, condition): the
                                     top-ranked sample + mean/sd over 5 samples
  eval/<ID>/per_residue_<cond>.csv   top-ranked sample: pLDDT, lDDT, Cα error
                                     after the TM-score superposition

Usage:
    python scripts/evaluate_structures.py
    python scripts/evaluate_structures.py --target T1183 --condition custom
"""
import argparse
import csv
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one

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
LDDT_RADIUS = 15.0
LDDT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)

_parser = MMCIFParser(QUIET=True)


def reference_path(casp_id: str) -> Path:
    pdb_id, chain = TARGET_REFS[casp_id]
    if chain:
        return RCSB_CHAIN_DIR / f"{pdb_id}_{chain}.cif"
    return RCSB_DIR / f"{pdb_id}.cif"


def eval_unit(casp_id: str):
    with open(ROOT / "data" / "metadata.csv") as fh:
        for row in csv.DictReader(fh):
            if row["casp_id"] == casp_id:
                lo, hi = row["eval_unit"].split("-")
                return int(lo), int(hi)
    raise KeyError(casp_id)


def ca_table(cif_path: Path):
    """resnum -> (one-letter aa, CA xyz, CA B-factor), first model, first chain.
    Keeps MSE (selenomethionine, HETATM in 8ORK) as M."""
    chain = next(iter(_parser.get_structure("s", str(cif_path))[0]))
    out = {}
    for res in chain:
        hetflag, resnum, _ = res.id
        if "CA" not in res or (hetflag.strip() and res.get_resname() != "MSE"):
            continue
        name = "MET" if res.get_resname() == "MSE" else res.get_resname()
        try:
            aa = index_to_one(three_to_index(name))
        except (KeyError, ValueError):
            aa = "X"
        out[resnum] = (aa, res["CA"].coord.astype(float), float(res["CA"].bfactor))
    return out


ONE_TO_THREE = {index_to_one(i): name for i, name in enumerate(
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}


def write_ca_pdb(path: Path, table, resnums):
    lines = [
        "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C"
        % (i + 1, ONE_TO_THREE.get(table[r][0], "UNK"), r, *table[r][1], table[r][2])
        for i, r in enumerate(resnums)
    ]
    path.write_text("\n".join(lines + ["TER", "END"]) + "\n")


def run_usalign(model, ref, resnums):
    """US-align on C-alpha-only PDBs restricted to `resnums`. Writing our own
    files (a) applies the CASP evaluation unit and (b) keeps MSE residues,
    which US-align silently drops when they are HETATM records (8ORK has 12)."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        write_ca_pdb(tmp / "model.pdb", model, resnums)
        write_ca_pdb(tmp / "ref.pdb", ref, resnums)
        matrix_file = tmp / "m.txt"
        out = subprocess.run(
            [str(USALIGN), str(tmp / "model.pdb"), str(tmp / "ref.pdb"),
             "-TMscore", "1", "-m", str(matrix_file)],
            capture_output=True, text=True, check=True,
        ).stdout
        n_ref = int(re.search(r"Length of Structure_2:\s*(\d+)", out).group(1))
        if n_ref != len(resnums):
            raise SystemExit(f"US-align read {n_ref} reference residues, expected {len(resnums)}")
        rows = [l.split() for l in matrix_file.read_text().splitlines() if re.match(r"^[012]\s", l)]
    t = np.array([float(r[1]) for r in rows])
    u = np.array([[float(x) for x in r[2:5]] for r in rows])
    tm = float(re.search(r"TM-score=\s*([\d.]+) \(normalized by length of Structure_2", out).group(1))
    rmsd = float(re.search(r"RMSD=\s*([\d.]+)", out).group(1))
    return tm, rmsd, u, t


def lddt_per_residue(model_xyz: np.ndarray, ref_xyz: np.ndarray) -> np.ndarray:
    d_ref = np.linalg.norm(ref_xyz[:, None] - ref_xyz[None], axis=-1)
    d_mod = np.linalg.norm(model_xyz[:, None] - model_xyz[None], axis=-1)
    mask = (d_ref < LDDT_RADIUS) & ~np.eye(len(ref_xyz), dtype=bool)
    diff = np.abs(d_ref - d_mod)
    preserved = np.mean([(diff < t) & mask for t in LDDT_THRESHOLDS], axis=0).sum(axis=1)
    n = mask.sum(axis=1)
    return np.where(n > 0, preserved / np.maximum(n, 1), np.nan)


def score_model(casp_id: str, model_path: Path, full_data_path: Path):
    ref_path = reference_path(casp_id)
    lo, hi = eval_unit(casp_id)
    ref, model = ca_table(ref_path), ca_table(model_path)
    resnums = [r for r in sorted(ref) if lo <= r <= hi and r in model]
    mismatches = [r for r in resnums if ref[r][0] != model[r][0]]
    if mismatches:
        raise SystemExit(f"[{casp_id}] residue identity mismatch model vs reference at {mismatches[:5]}")

    tm, rmsd, u, t = run_usalign(model, ref, resnums)
    ref_xyz = np.array([ref[r][1] for r in resnums])
    mod_xyz = np.array([model[r][1] for r in resnums])
    lddt = lddt_per_residue(mod_xyz, ref_xyz)
    ca_err = np.linalg.norm((mod_xyz @ u.T + t) - ref_xyz, axis=1)

    plddt_all = np.array([model[r][2] for r in sorted(model)])
    pae = np.array(json.loads(full_data_path.read_text())["pae"], dtype=float)

    per_res = [
        {"resnum": r, "aa": model[r][0], "plddt": round(model[r][2], 2),
         "lddt": round(float(l), 4), "ca_error": round(float(e), 3)}
        for r, l, e in zip(resnums, lddt, ca_err)
    ]
    # residues outside the evaluation unit / unresolved in the reference
    # still get their pLDDT, with blank accuracy columns
    scored = set(resnums)
    per_res += [
        {"resnum": r, "aa": model[r][0], "plddt": round(model[r][2], 2), "lddt": "", "ca_error": ""}
        for r in model if r not in scored
    ]
    per_res.sort(key=lambda d: d["resnum"])

    return {
        "tm_score": tm,
        "rmsd": rmsd,
        "lddt": float(np.nanmean(lddt)),
        "mean_plddt": float(plddt_all.mean()),
        "mean_plddt_eval_unit": float(np.mean([model[r][2] for r in resnums])),
        "mean_pae": float(pae.mean()),
        "n_scored_residues": len(resnums),
    }, per_res


def evaluate_condition(casp_id: str, condition: str):
    cond_dir = ROOT / f"af3_{condition}" / casp_id
    meta_path = cond_dir / "job_metadata.json"
    out_dir = cond_dir / "server_output"
    if not meta_path.exists() or not out_dir.exists():
        print(f"[{casp_id}/{condition}] no imported AF3 output yet, skipping "
              "(run scripts/import_af3_results.py)")
        return [], None
    meta = json.loads(meta_path.read_text())

    rows, top_per_res = [], None
    for model_path in sorted(out_dir.glob("*_model_*.cif")):
        idx = int(model_path.stem.rsplit("_", 1)[1])
        stem = model_path.stem.rsplit("_model_", 1)[0]
        summary = json.loads((out_dir / f"{stem}_summary_confidences_{idx}.json").read_text())
        scores, per_res = score_model(casp_id, model_path, out_dir / f"{stem}_full_data_{idx}.json")
        rows.append({
            "casp_id": casp_id, "condition": condition, "sample": idx,
            "top_ranked": idx == meta["top_ranked_sample"],
            "ranking_score": summary["ranking_score"], "ptm": summary["ptm"],
            **scores,
        })
        if idx == meta["top_ranked_sample"]:
            top_per_res = per_res
    return rows, top_per_res


def write_csv(path: Path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows):
    top = next(r for r in rows if r["top_ranked"])
    out = {k: top[k] for k in ("casp_id", "condition", "tm_score", "rmsd", "lddt",
                               "mean_plddt", "mean_pae", "ptm", "ranking_score", "n_scored_residues")}
    for metric in ("tm_score", "lddt", "rmsd"):
        vals = np.array([r[metric] for r in rows])
        out[f"{metric}_mean5"] = round(float(vals.mean()), 4)
        out[f"{metric}_sd5"] = round(float(vals.std(ddof=1)), 4) if len(vals) > 1 else 0.0
        out[f"{metric}_best5"] = round(float(vals.min() if metric == "rmsd" else vals.max()), 4)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=list(TARGET_REFS))
    parser.add_argument("--condition", choices=CONDITIONS)
    args = parser.parse_args()

    if not USALIGN.exists():
        raise SystemExit(f"{USALIGN} not found -- run scripts/build_usalign.sh first")

    targets = [args.target] if args.target else list(TARGET_REFS)
    conditions = [args.condition] if args.condition else CONDITIONS

    all_rows, summaries = [], []
    for casp_id in targets:
        for condition in conditions:
            rows, per_res = evaluate_condition(casp_id, condition)
            if not rows:
                continue
            all_rows += rows
            summaries.append(summarize(rows))
            write_csv(EVAL_DIR / casp_id / f"per_residue_{condition}.csv", per_res)
            s = summaries[-1]
            print(f"[{casp_id}/{condition}] top-ranked: TM={s['tm_score']:.3f} RMSD={s['rmsd']:.2f} "
                  f"lDDT={s['lddt']:.3f} pLDDT={s['mean_plddt']:.1f} | "
                  f"5-sample TM {s['tm_score_mean5']:.3f}±{s['tm_score_sd5']:.3f}")

    if not all_rows:
        print("no AF3 models found -- nothing to score yet")
        return

    for r in all_rows:
        for k, v in r.items():
            if isinstance(v, float):
                r[k] = round(v, 4)
    for s in summaries:
        for k, v in s.items():
            if isinstance(v, float):
                s[k] = round(v, 4)

    # merge with rows for targets/conditions not re-scored in this call
    def merge(path, new, key):
        keep = []
        if path.exists():
            with open(path) as fh:
                done = {key(r) for r in new}
                keep = [r for r in csv.DictReader(fh) if key(r) not in done]
        merged = sorted(keep + new, key=lambda r: (r["casp_id"], r["condition"], int(r.get("sample", 0))))
        write_csv(path, merged, list(new[0]))
        print(f"wrote {path.relative_to(ROOT)} ({len(merged)} rows)")

    merge(EVAL_DIR / "scores_all_samples.csv", all_rows, lambda r: (r["casp_id"], r["condition"]))
    merge(EVAL_DIR / "scores_summary.csv", summaries, lambda r: (r["casp_id"], r["condition"]))


if __name__ == "__main__":
    main()
