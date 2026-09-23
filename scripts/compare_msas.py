#!/usr/bin/env python3
"""Compare our custom MSA with the MSA the AF3 server actually built for the
default condition (both are returned in the AF3 download, see
af3_default/<ID>/server_output/msas/).

Per target and MSA source, computes depth, Neff/L (80% identity clustering,
same definition as scripts/msa_qc.py), mean hit coverage and identity, and
per-column coverage (how many non-query sequences have a residue at each
query position). AF3 feeds a monomer both its unpaired MSA and its
(UniProt) "paired" MSA, so the default condition is also reported as the
union of the two ("af3_default_all") -- that is what the model actually saw.

Writes:
  eval/msa_comparison.csv                  one row per (target, msa source)
  eval/<ID>/msa_column_coverage.csv        per-position coverage, both sources

Usage:
    python scripts/compare_msas.py
"""
import csv
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from msa_qc import compute_neff, match_columns, pairwise_identity, parse_a3m  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = ROOT / "eval"
TARGETS = ["T1183", "T1112", "T1122"]


def msa_stats(path: Path):
    cols = [match_columns(s) for _, s in parse_a3m(path)]
    L = len(cols[0])
    if any(len(c) != L for c in cols):
        raise SystemExit(f"{path}: rows have unequal match-column length")
    hits = cols[1:]
    col_cov = np.array([sum(1 for h in hits if h[i] != "-") for i in range(L)])
    neff = compute_neff(cols, 0.8, chunk_size=100)
    return {
        "depth": len(cols),
        "neff": round(neff, 2),
        "neff_per_l": round(neff / L, 4),
        "mean_coverage": round(float(np.mean([sum(c != "-" for c in h) / L for h in hits])), 4) if hits else 0.0,
        "mean_identity": round(float(np.mean([pairwise_identity(cols[0], h) for h in hits])), 4) if hits else 0.0,
        "min_column_coverage": int(col_cov.min()),
    }, col_cov


def concat_a3m(paths, out: Path):
    """Union of several A3Ms that share a query: keep one query row, then all hits."""
    records = []
    for i, p in enumerate(paths):
        recs = list(parse_a3m(p))
        records += recs if i == 0 else recs[1:]
    out.write_text("".join(f">{h}\n{s}\n" for h, s in records))
    return out


def main():
    tmp = Path(tempfile.mkdtemp())
    rows = []
    for casp_id in TARGETS:
        server_msas = ROOT / "af3_default" / casp_id / "server_output" / "msas"
        sources = {
            "custom": ROOT / "msa" / casp_id / f"{casp_id}_custom.a3m",
            "af3_default_unpaired": next(server_msas.glob("*_unpaired_msa_*.a3m"), None),
        }
        if sources["af3_default_unpaired"] is None:
            print(f"[{casp_id}] no AF3 default MSA imported yet, skipping")
            continue
        paired = next(server_msas.glob("*_paired_msa_*.a3m"), None)
        n_paired = sum(1 for _ in parse_a3m(paired)) if paired else 0
        if paired:
            sources["af3_default_all"] = concat_a3m(
                [sources["af3_default_unpaired"], paired], tmp / f"{casp_id}_af3_all.a3m")

        coverage = {}
        for source, path in sources.items():
            stats, coverage[source] = msa_stats(path)
            rows.append({"casp_id": casp_id, "msa_source": source, **stats,
                         "af3_paired_depth": n_paired if source != "custom" else "",
                         "path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT)
                         else "unpaired + paired server MSAs"})
            print(f"[{casp_id}/{source}] depth={stats['depth']} Neff/L={stats['neff_per_l']}")

        out = EVAL_DIR / casp_id / "msa_column_coverage.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["resnum", *coverage])
            for i in range(len(coverage["custom"])):
                w.writerow([i + 1, *(int(c[i]) for c in coverage.values())])

    out = EVAL_DIR / "msa_comparison.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
