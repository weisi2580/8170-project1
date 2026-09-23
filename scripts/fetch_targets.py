#!/usr/bin/env python3
"""Fetch/verify CASP15 target FASTAs + RCSB reference structures, and extract the
evaluation chain from each multi-chain reference.

This script documents (and makes reproducible) the download that was originally
done by hand through the CASP15 and RCSB PDB websites. Re-running it with
--download will re-fetch the FASTA/CIF files from RCSB; without --download it
only (re)builds the chain-filtered structures and data/metadata.csv from
whatever is already in data/CASP15 and data/RCSB.

Usage:
    python scripts/fetch_targets.py                # build chain-filtered PDBs + metadata.csv
    python scripts/fetch_targets.py --download      # also (re)download FASTA/CIF from RCSB first
"""
import argparse
import csv
import datetime as dt
from pathlib import Path

import requests
from Bio.PDB import MMCIFParser, MMCIFIO, Select

ROOT = Path(__file__).resolve().parents[1]
CASP_DIR = ROOT / "data" / "CASP15"
RCSB_DIR = ROOT / "data" / "RCSB"
CHAIN_DIR = ROOT / "data" / "RCSB_chain"
METADATA_CSV = ROOT / "data" / "metadata.csv"

# (CASP target ID, PDB ID, chain to keep, CASP15 difficulty class, evaluation unit)
# Difficulty class and evaluation-unit residue range are copied from the
# official CASP15 domain table (checked 2026-09-23):
# https://predictioncenter.org/casp15/domains_summary.cgi
TARGETS = [
    {"casp_id": "T1183", "pdb_id": "8IFX", "chain": "B", "difficulty": "TBM-easy", "eval_unit": "1-195"},
    {"casp_id": "T1112", "pdb_id": "8ORK", "chain": None, "difficulty": "FM/TBM", "eval_unit": "1-460"},  # single-chain
    {"casp_id": "T1122", "pdb_id": "8BBT", "chain": "A", "difficulty": "FM", "eval_unit": "4-237"},
]

RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"


class ChainSelector(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def download_cif(pdb_id: str) -> Path:
    dest = RCSB_DIR / f"{pdb_id}.cif"
    RCSB_DIR.mkdir(parents=True, exist_ok=True)
    resp = requests.get(RCSB_CIF_URL.format(pdb_id=pdb_id), timeout=30)
    resp.raise_for_status()
    dest.write_text(resp.text)
    return dest


def extract_chain(pdb_id: str, chain_id: str) -> Path:
    """Write a single-chain copy of an mmCIF file to data/RCSB_chain/."""
    CHAIN_DIR.mkdir(parents=True, exist_ok=True)
    src = RCSB_DIR / f"{pdb_id}.cif"
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, src)

    io = MMCIFIO()
    io.set_structure(structure)
    dest = CHAIN_DIR / f"{pdb_id}_{chain_id}.cif"
    io.save(str(dest), ChainSelector(chain_id))
    return dest


def read_fasta_length(fasta_path: Path) -> int:
    seq = ""
    for line in fasta_path.read_text().splitlines():
        if not line.startswith(">"):
            seq += line.strip()
    return len(seq)


def build_metadata(rows):
    METADATA_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(METADATA_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["casp_id", "length", "pdb_id", "chain", "download_date", "casp_difficulty_class", "eval_unit"],
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="re-download CIF files from RCSB before processing")
    args = parser.parse_args()

    today = dt.date.today().isoformat()
    rows = []

    for target in TARGETS:
        casp_id, pdb_id, chain, difficulty = (
            target["casp_id"],
            target["pdb_id"],
            target["chain"],
            target["difficulty"],
        )

        if args.download:
            print(f"downloading {pdb_id}.cif ...")
            download_cif(pdb_id)

        fasta_path = CASP_DIR / f"{casp_id}.fasta"
        length = read_fasta_length(fasta_path)

        if chain is not None:
            out = extract_chain(pdb_id, chain)
            print(f"{pdb_id}: wrote single-chain structure {out.relative_to(ROOT)} (chain {chain})")
        else:
            print(f"{pdb_id}: single-chain entry, no extraction needed")

        rows.append(
            {
                "casp_id": casp_id,
                "length": length,
                "pdb_id": pdb_id,
                "chain": chain or "(single chain)",
                "download_date": today,
                "casp_difficulty_class": difficulty,
                "eval_unit": target["eval_unit"],
            }
        )

    build_metadata(rows)
    print(f"wrote {METADATA_CSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
