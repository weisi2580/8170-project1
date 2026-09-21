#!/usr/bin/env python3
"""MSA quality-control metrics and figures for a filtered custom A3M.

For the target's msa/<ID>/<ID>_custom.a3m, computes:
  1. depth and Neff/L (effective sequence count at 80% identity, normalized
     by query length)
  2. per-hit query coverage distribution
  3. per-column gap fraction
  4. query-identity distribution across hits
  5. per-column Shannon entropy (conservation) + a sequence logo

Writes figures/<ID>/msa_qc_*.png and appends one row to
eval/msa_qc_summary.csv (one row per target, so all 3 targets end up
compared in a single table).

Usage:
    python scripts/msa_qc.py --target T1183
    python scripts/msa_qc.py --target T1183 --neff-identity 0.8
"""
import argparse
import csv
import math
from collections import Counter
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
MSA_DIR = ROOT / "msa"
FIG_DIR = ROOT / "figures"
EVAL_DIR = ROOT / "eval"

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

sns.set_theme(style="whitegrid")


def parse_a3m(path: Path):
    header, seq_lines = None, []
    for line in path.read_text().splitlines():
        line = line.replace("\x00", "")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_lines)
            header, seq_lines = line[1:], []
        else:
            seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def match_columns(seq: str) -> str:
    return "".join(c for c in seq if c.isupper() or c == "-")


def pairwise_identity(a: str, b: str) -> float:
    """Fraction identical over columns where BOTH sequences have a residue
    (same convention used in scripts/filter_msa.py)."""
    n_both = sum(1 for x, y in zip(a, b) if x != "-" and y != "-")
    n_match = sum(1 for x, y in zip(a, b) if x != "-" and y != "-" and x == y)
    return n_match / n_both if n_both else 0.0


def compute_neff(cols, identity_threshold: float, chunk_size: int = 500) -> float:
    """Standard clustering Neff: each sequence's weight is 1 / (# sequences
    within identity_threshold of it, including itself). Vectorized with numpy
    and row-chunked so memory stays bounded for large alignments."""
    n = len(cols)
    L = len(cols[0])
    # encode: 0 = gap, 1.. = ord(char); any two non-gap equal chars match
    enc = np.array([[ord(c) if c != "-" else -1 for c in row] for row in cols], dtype=np.int16)
    is_res = enc != -1

    cluster_sizes = np.zeros(n, dtype=np.float64)
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        chunk_enc = enc[start:end]  # (c, L)
        chunk_res = is_res[start:end]  # (c, L)
        # both sequences have a residue at each position, for normalizing identity
        both = chunk_res[:, None, :] & is_res[None, :, :]  # (c, n, L)
        match = (chunk_enc[:, None, :] == enc[None, :, :]) & chunk_res[:, None, :]  # (c, n, L)
        n_match = match.sum(axis=2).astype(np.float64)
        n_both = both.sum(axis=2).astype(np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            ident = np.where(n_both > 0, n_match / n_both, 0.0)
        cluster_sizes[start:end] = (ident >= identity_threshold).sum(axis=1)

    weights = 1.0 / cluster_sizes
    return float(weights.sum())


def column_entropy(columns_by_position):
    entropies = []
    for col in columns_by_position:
        counts = Counter(c for c in col if c != "-" and c in AMINO_ACIDS)
        total = sum(counts.values())
        if total == 0:
            entropies.append(0.0)
            continue
        ent = -sum((c / total) * math.log2(c / total) for c in counts.values())
        entropies.append(ent)
    return entropies


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument("--neff-identity", type=float, default=0.8)
    args = parser.parse_args()
    target = args.target

    a3m_path = MSA_DIR / target / f"{target}_custom.a3m"
    fig_dir = FIG_DIR / target
    fig_dir.mkdir(parents=True, exist_ok=True)
    EVAL_DIR.mkdir(parents=True, exist_ok=True)

    records = list(parse_a3m(a3m_path))
    headers = [h for h, _ in records]
    cols = [match_columns(s) for _, s in records]
    query_cols = cols[0]
    L = len(query_cols)
    depth = len(records)

    # 1. depth & Neff/L
    neff = compute_neff(cols, args.neff_identity)
    neff_per_l = neff / L

    # 2. per-hit coverage distribution (hits only, excludes query row 0)
    coverages = [
        sum(1 for c in hit if c != "-") / L for hit in cols[1:]
    ]
    plt.figure(figsize=(6, 4))
    sns.histplot(coverages, bins=20)
    plt.xlabel("Query coverage")
    plt.ylabel("Number of hits")
    plt.title(f"{target}: MSA hit coverage distribution")
    plt.tight_layout()
    plt.savefig(fig_dir / "msa_qc_coverage.png", dpi=150)
    plt.close()

    # 3. per-column gap fraction
    gap_fraction = [
        sum(1 for hit in cols if hit[i] == "-") / depth for i in range(L)
    ]
    plt.figure(figsize=(10, 3))
    plt.plot(range(1, L + 1), gap_fraction)
    plt.xlabel("Query column")
    plt.ylabel("Gap fraction")
    plt.title(f"{target}: per-column gap fraction ({depth} sequences)")
    plt.tight_layout()
    plt.savefig(fig_dir / "msa_qc_gap_fraction.png", dpi=150)
    plt.close()

    # 4. query-identity distribution
    identities = [pairwise_identity(query_cols, hit) for hit in cols[1:]]
    plt.figure(figsize=(6, 4))
    sns.histplot(identities, bins=20)
    plt.xlabel("Identity to query")
    plt.ylabel("Number of hits")
    plt.title(f"{target}: hit identity-to-query distribution")
    plt.tight_layout()
    plt.savefig(fig_dir / "msa_qc_identity.png", dpi=150)
    plt.close()

    # 5. conservation / entropy + sequence logo
    columns_by_position = [[hit[i] for hit in cols] for i in range(L)]
    entropies = column_entropy(columns_by_position)
    plt.figure(figsize=(10, 3))
    plt.plot(range(1, L + 1), entropies)
    plt.xlabel("Query column")
    plt.ylabel("Shannon entropy (bits)")
    plt.title(f"{target}: per-column conservation (entropy)")
    plt.tight_layout()
    plt.savefig(fig_dir / "msa_qc_entropy.png", dpi=150)
    plt.close()

    try:
        import logomaker
        import pandas as pd

        freq_rows = []
        for col in columns_by_position:
            counts = Counter(c for c in col if c != "-" and c in AMINO_ACIDS)
            total = sum(counts.values()) or 1
            freq_rows.append({aa: counts.get(aa, 0) / total for aa in AMINO_ACIDS})
        freq_df = pd.DataFrame(freq_rows)
        plt.figure(figsize=(max(10, L * 0.08), 3))
        logomaker.Logo(freq_df, figsize=(max(10, L * 0.08), 3))
        plt.title(f"{target}: sequence logo (custom MSA)")
        plt.tight_layout()
        plt.savefig(fig_dir / "msa_qc_logo.png", dpi=150)
        plt.close()
    except Exception as e:  # pragma: no cover - logo is a nice-to-have
        print(f"[{target}] skipped sequence logo: {e}")

    # summary CSV (append/replace row for this target)
    summary_path = EVAL_DIR / "msa_qc_summary.csv"
    fieldnames = ["casp_id", "length", "depth", "neff", "neff_per_l", "mean_coverage", "mean_identity"]
    row = {
        "casp_id": target,
        "length": L,
        "depth": depth,
        "neff": round(neff, 2),
        "neff_per_l": round(neff_per_l, 4),
        "mean_coverage": round(sum(coverages) / len(coverages), 4) if coverages else 0.0,
        "mean_identity": round(sum(identities) / len(identities), 4) if identities else 0.0,
    }
    existing_rows = []
    if summary_path.exists():
        with open(summary_path) as fh:
            existing_rows = [r for r in csv.DictReader(fh) if r["casp_id"] != target]
    existing_rows.append(row)
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(existing_rows)

    print(f"[{target}] depth={depth} Neff={neff:.1f} Neff/L={neff_per_l:.3f}")
    print(f"[{target}] figures written to {fig_dir.relative_to(ROOT)}")
    print(f"[{target}] summary row written to {summary_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
