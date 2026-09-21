#!/usr/bin/env python3
"""Standardize the AF3 job metadata schema and generate the submission checklist.

The public AlphaFold3 server (https://alphafoldserver.com) has no batch API for
this use case, so the 6 jobs (3 targets x {custom, default}) must be submitted
by hand in the browser. This script does the parts that ARE scriptable:

  --prep    write the exact per-target input (sequence, and for the custom
            condition the path to the filtered A3M) plus a checklist to
            af3_custom/<ID>/ and af3_default/<ID>/, so nothing has to be
            retyped in the browser.
  --log     after a job is submitted, record its metadata (job id, submission
            timestamp, settings) into af3_custom/<ID>/job_metadata.json or
            af3_default/<ID>/job_metadata.json with a fixed schema, so all 6
            runs are comparable.

Usage:
    python scripts/log_af3_job.py --prep
    python scripts/log_af3_job.py --log --target T1183 --condition custom \
        --job-id abc123 --submitted-at 2026-09-22T10:00:00 --model-version AF3-server-2024.5
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASP_DIR = ROOT / "data" / "CASP15"
MSA_DIR = ROOT / "msa"

TARGETS = ["T1183", "T1112", "T1122"]
CONDITIONS = ["custom", "default"]

SCHEMA_FIELDS = [
    "casp_id",
    "condition",  # "custom" | "default"
    "job_id",
    "submitted_at",  # ISO 8601
    "model_version",
    "msa_input",  # path to the a3m used, or "AF3 server default MSA"
    "notes",
]


def read_seq(casp_id: str) -> str:
    fasta = CASP_DIR / f"{casp_id}.fasta"
    return "".join(
        l.strip() for l in fasta.read_text().splitlines() if not l.startswith(">")
    )


def prep():
    for casp_id in TARGETS:
        seq = read_seq(casp_id)

        custom_dir = ROOT / "af3_custom" / casp_id
        custom_dir.mkdir(parents=True, exist_ok=True)
        a3m_path = MSA_DIR / casp_id / f"{casp_id}_custom.a3m"
        (custom_dir / "input_sequence.fasta").write_text(f">{casp_id}\n{seq}\n")
        checklist = f"""AF3 submission checklist -- {casp_id} (custom MSA condition)
1. Go to https://alphafoldserver.com and start a new job.
2. Paste the sequence from input_sequence.fasta (or upload it).
3. Under MSA/templates, upload the custom alignment:
     {a3m_path.relative_to(ROOT)}
   (must exist -- run scripts/build_msa.py + scripts/filter_msa.py first)
4. Leave every other setting at the same value used for the default-MSA run.
5. Submit. Record the job ID + timestamp with:
     python scripts/log_af3_job.py --log --target {casp_id} --condition custom \\
         --job-id <ID> --submitted-at <ISO8601> --model-version <version>
6. Once complete, download the top-ranked model (.cif/.pdb) and confidence
   JSON into this folder.
"""
        (custom_dir / "SUBMIT_CHECKLIST.txt").write_text(checklist)

        default_dir = ROOT / "af3_default" / casp_id
        default_dir.mkdir(parents=True, exist_ok=True)
        (default_dir / "input_sequence.fasta").write_text(f">{casp_id}\n{seq}\n")
        checklist_default = f"""AF3 submission checklist -- {casp_id} (default MSA condition)
1. Go to https://alphafoldserver.com and start a new job.
2. Paste the raw sequence from input_sequence.fasta -- do NOT upload a custom MSA,
   let AF3 build its own alignment.
3. Leave every other setting identical to the custom-MSA run for {casp_id}.
4. Submit. Record the job ID + timestamp with:
     python scripts/log_af3_job.py --log --target {casp_id} --condition default \\
         --job-id <ID> --submitted-at <ISO8601> --model-version <version>
5. Once complete, download the top-ranked model (.cif/.pdb) and confidence
   JSON into this folder.
"""
        (default_dir / "SUBMIT_CHECKLIST.txt").write_text(checklist_default)

        print(f"[{casp_id}] wrote input + checklist for both conditions")


def log(args):
    condition_dir = ROOT / f"af3_{args.condition}" / args.target
    condition_dir.mkdir(parents=True, exist_ok=True)

    msa_input = (
        str((MSA_DIR / args.target / f"{args.target}_custom.a3m").relative_to(ROOT))
        if args.condition == "custom"
        else "AF3 server default MSA"
    )

    record = {
        "casp_id": args.target,
        "condition": args.condition,
        "job_id": args.job_id,
        "submitted_at": args.submitted_at,
        "model_version": args.model_version,
        "msa_input": msa_input,
        "notes": args.notes or "",
    }
    out_path = condition_dir / "job_metadata.json"
    out_path.write_text(json.dumps(record, indent=2) + "\n")
    print(f"wrote {out_path.relative_to(ROOT)}")


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prep", action="store_true")
    mode.add_argument("--log", action="store_true")
    parser.add_argument("--target", choices=TARGETS)
    parser.add_argument("--condition", choices=CONDITIONS)
    parser.add_argument("--job-id")
    parser.add_argument("--submitted-at", help="ISO 8601 timestamp")
    parser.add_argument("--model-version")
    parser.add_argument("--notes", default="")
    args = parser.parse_args()

    if args.prep:
        prep()
    else:
        missing = [
            name
            for name, val in [
                ("--target", args.target),
                ("--condition", args.condition),
                ("--job-id", args.job_id),
                ("--submitted-at", args.submitted_at),
                ("--model-version", args.model_version),
            ]
            if not val
        ]
        if missing:
            parser.error(f"--log requires {', '.join(missing)}")
        log(args)


if __name__ == "__main__":
    main()
