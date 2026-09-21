#!/usr/bin/env python3
"""Filter a raw A3M (from scripts/build_msa.py) into the final custom MSA.

Filters applied, in order (each threshold is a CLI flag so choices can be
justified/tuned in the report):
  1. query coverage  -- drop hits aligning to fewer than --min-coverage of
                         query (non-gap) columns.
  2. redundancy       -- greedily drop hits >= --max-identity identical to a
                         hit already kept (simple single-linkage redundancy
                         filter, similar in spirit to `hhfilter -id`).
  3. min identity      -- drop hits below --min-identity to the query
                         (removes likely spurious/very-remote hits).

Also validates the output: row 0 must equal the query FASTA exactly, every
row must be the same alignment length, and only A3M-legal characters
(upper/lower letters, '-', '.') are present.

Usage:
    python scripts/filter_msa.py --target T1183
    python scripts/filter_msa.py --target T1183 --min-coverage 0.5 --max-identity 0.9 --min-identity 0.2
"""
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASP_DIR = ROOT / "data" / "CASP15"
MSA_DIR = ROOT / "msa"

A3M_LEGAL = set("ACDEFGHIKLMNPQRSTVWYXBZJUOacdefghiklmnpqrstvwyxbzjuo-.*\n")


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
    """Return only the match-state (query-aligned) columns of an A3M row:
    uppercase letters and '-', dropping lowercase insertions."""
    return "".join(c for c in seq if c.isupper() or c == "-")


def coverage(query_cols: str, hit_cols: str) -> float:
    aligned = sum(1 for q, h in zip(query_cols, hit_cols) if h != "-")
    return aligned / len(query_cols)


def identity(query_cols: str, hit_cols: str) -> float:
    """Fraction identical over columns where BOTH sequences have a residue
    (same convention as pairwise_identity below, and as scripts/msa_qc.py)."""
    return pairwise_identity(query_cols, hit_cols)


def pairwise_identity(a: str, b: str) -> float:
    n_both = sum(1 for x, y in zip(a, b) if x != "-" and y != "-")
    n_match = sum(1 for x, y in zip(a, b) if x != "-" and y != "-" and x == y)
    return n_match / n_both if n_both else 0.0


def filter_msa(records, min_coverage, max_identity, min_identity):
    query_header, query_seq = records[0]
    query_cols = match_columns(query_seq)

    kept = [(query_header, query_seq, query_cols)]
    for header, seq in records[1:]:
        hit_cols = match_columns(seq)
        if len(hit_cols) != len(query_cols):
            continue  # malformed row, skip
        cov = coverage(query_cols, hit_cols)
        if cov < min_coverage:
            continue
        idty = identity(query_cols, hit_cols)
        if idty < min_identity:
            continue
        # redundancy filter against everything kept so far
        if any(pairwise_identity(hit_cols, k[2]) >= max_identity for k in kept[1:]):
            continue
        kept.append((header, seq, hit_cols))

    return kept


def validate(kept, query_fasta_seq: str):
    assert kept, "no rows survived filtering"
    header0, seq0, cols0 = kept[0]
    assert match_columns(seq0).replace("-", "") == query_fasta_seq, (
        "row 0 does not match the query FASTA exactly"
    )
    lengths = {len(match_columns(seq)) for _, seq, _ in kept}
    assert len(lengths) == 1, f"inconsistent alignment lengths: {lengths}"
    for _, seq, _ in kept:
        assert set(seq) <= A3M_LEGAL, "illegal character in A3M row"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument("--min-coverage", type=float, default=0.5)
    parser.add_argument("--max-identity", type=float, default=0.9)
    parser.add_argument("--min-identity", type=float, default=0.2)
    args = parser.parse_args()

    raw_path = MSA_DIR / args.target / f"{args.target}_raw.a3m"
    out_path = MSA_DIR / args.target / f"{args.target}_custom.a3m"
    fasta_path = CASP_DIR / f"{args.target}.fasta"

    query_fasta_seq = "".join(
        l.strip() for l in fasta_path.read_text().splitlines() if not l.startswith(">")
    )

    records = list(parse_a3m(raw_path))
    kept = filter_msa(records, args.min_coverage, args.max_identity, args.min_identity)
    validate(kept, query_fasta_seq)

    with open(out_path, "w") as fh:
        for header, seq, _ in kept:
            fh.write(f">{header}\n{seq}\n")

    print(
        f"[{args.target}] {len(records)} raw rows -> {len(kept)} kept "
        f"(min_coverage={args.min_coverage}, max_identity={args.max_identity}, "
        f"min_identity={args.min_identity})"
    )
    print(f"[{args.target}] wrote {out_path.relative_to(ROOT)} (validated)")


if __name__ == "__main__":
    main()
