#!/usr/bin/env python3
"""Build a custom MSA per target using the public ColabFold MMseqs2 search server.

This project doesn't have local access to UniRef30 / BFD-lite / ColabFold-envdb
(hundreds of GB) or a local mmseqs2 binary, so instead of running
`colabfold_search` against local databases, this script talks to the same
remote MMseqs2 search API that the ColabFold notebook (api.colabfold.com) uses.
It is the "accessible via ColabFold's remote MSA server" path mentioned in
project_plan.md Step 3. This hits a shared public, rate-limited, non-commercial
research server -- be polite (this script already paces itself), and don't run
it in parallel for multiple targets.

For each target it:
  1. submits the FASTA sequence in `mode=env` (paired UniRef100 search +
     environmental/metagenomic search, matching "UniRef30 + BFD-lite/MGnify"),
  2. polls until the job is COMPLETE,
  3. downloads and unpacks the result tarball, saving the raw per-database
     A3M files under msa/<ID>/raw/,
  4. concatenates uniref.a3m + bfd.mgnify30.metaeuk30.smag30.a3m into a single
     msa/<ID>/<ID>_raw.a3m (query row first, de-duplicated by sequence) that
     scripts/filter_msa.py then filters into the final custom MSA.

Usage:
    python scripts/build_msa.py                 # all 3 targets
    python scripts/build_msa.py --target T1183   # just one
"""
import argparse
import io
import tarfile
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
CASP_DIR = ROOT / "data" / "CASP15"
MSA_DIR = ROOT / "msa"

API_BASE = "https://api.colabfold.com"
TARGETS = ["T1183", "T1112", "T1122"]

POLL_INTERVAL_SEC = 10
MAX_POLL_SEC = 20 * 60


def submit_job(fasta_text: str) -> str:
    resp = requests.post(
        f"{API_BASE}/ticket/msa",
        data={"q": fasta_text, "mode": "env"},
        timeout=30,
    )
    resp.raise_for_status()
    payload = resp.json()
    if "id" not in payload:
        raise RuntimeError(f"unexpected submit response: {payload}")
    return payload["id"]


def wait_for_job(job_id: str) -> None:
    waited = 0
    while waited < MAX_POLL_SEC:
        resp = requests.get(f"{API_BASE}/ticket/{job_id}", timeout=15)
        resp.raise_for_status()
        status = resp.json().get("status")
        if status == "COMPLETE":
            return
        if status == "ERROR":
            raise RuntimeError(f"ColabFold server reported ERROR for job {job_id}")
        time.sleep(POLL_INTERVAL_SEC)
        waited += POLL_INTERVAL_SEC
    raise TimeoutError(f"job {job_id} did not complete within {MAX_POLL_SEC}s")


def download_result(job_id: str, dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    resp = requests.get(f"{API_BASE}/result/download/{job_id}", timeout=60)
    resp.raise_for_status()
    with tarfile.open(fileobj=io.BytesIO(resp.content), mode="r:gz") as tar:
        tar.extractall(dest_dir)


def parse_a3m_records(path: Path):
    """Yield (header, sequence) pairs from an A3M file. The ColabFold server
    occasionally pads the last sequence line with a stray NUL byte; strip it
    (and any other control characters) here so downstream files are clean."""
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


def merge_a3m(raw_dir: Path, casp_id: str, out_path: Path) -> int:
    """Concatenate the UniRef + environmental A3M hits into one raw A3M,
    query first, de-duplicated by aligned sequence. Returns number of rows written."""
    sources = [
        raw_dir / "uniref.a3m",
        raw_dir / "bfd.mgnify30.metaeuk30.smag30.a3m",
    ]
    seen_seqs = set()
    rows = []
    query_row = None

    for src in sources:
        if not src.exists():
            continue
        for header, seq in parse_a3m_records(src):
            if query_row is None and not rows:
                # first record of the first available source is the query itself
                query_row = (casp_id, seq)
                seen_seqs.add(seq)
                continue
            if seq in seen_seqs:
                continue
            seen_seqs.add(seq)
            rows.append((header, seq))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        if query_row:
            fh.write(f">{query_row[0]}\n{query_row[1]}\n")
        for header, seq in rows:
            fh.write(f">{header}\n{seq}\n")
    return len(rows) + (1 if query_row else 0)


def build_one(casp_id: str) -> None:
    fasta_path = CASP_DIR / f"{casp_id}.fasta"
    fasta_lines = fasta_path.read_text().splitlines()
    seq = "".join(l.strip() for l in fasta_lines if not l.startswith(">"))
    fasta_text = f">{casp_id}\n{seq}\n"

    print(f"[{casp_id}] submitting job ({len(seq)} aa) ...")
    job_id = submit_job(fasta_text)
    print(f"[{casp_id}] job id {job_id}, polling ...")
    wait_for_job(job_id)

    raw_dir = MSA_DIR / casp_id / "raw"
    print(f"[{casp_id}] downloading result ...")
    download_result(job_id, raw_dir)

    out_path = MSA_DIR / casp_id / f"{casp_id}_raw.a3m"
    n = merge_a3m(raw_dir, casp_id, out_path)
    print(f"[{casp_id}] wrote {out_path.relative_to(ROOT)} ({n} rows incl. query)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=TARGETS, help="only build for this target")
    args = parser.parse_args()

    targets = [args.target] if args.target else TARGETS
    for i, casp_id in enumerate(targets):
        build_one(casp_id)
        if i < len(targets) - 1:
            time.sleep(15)  # be polite to the shared public server between jobs


if __name__ == "__main__":
    main()
