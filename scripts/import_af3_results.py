#!/usr/bin/env python3
"""Import the AlphaFold Server download zips (Step 5 output) into the repo layout.

The AF3 server gives one zip per job, named fold_<target>_<condition>_msa.zip.
Each holds 5 ranked samples (model_0..4), their summary + full confidence
JSONs, the MSA(s) the server actually used, the template hits, and the job
request. This script:

  1. unpacks af3_raw/fold_t<ID>_<cond>_msa.zip -> af3_<cond>/<ID>/server_output/
  2. checks the job request's sequence is exactly data/CASP15/<ID>.fasta,
     and (custom condition) that the MSA the server used is byte-for-byte
     msa/<ID>/<ID>_custom.a3m
  3. copies the top-ranked sample (highest ranking_score) to
       af3_<cond>/<ID>/<ID>_model.cif
       af3_<cond>/<ID>/<ID>_confidence.json          (full_data: pLDDT, PAE)
       af3_<cond>/<ID>/<ID>_summary_confidences.json (pTM, ranking_score, ...)
  4. writes af3_<cond>/<ID>/job_metadata.json (same schema as
     scripts/log_af3_job.py --log, plus seed/template/sample fields)

Usage:
    python scripts/import_af3_results.py            # all zips in af3_raw/
    python scripts/import_af3_results.py --zip-dir some/other/dir
"""
import argparse
import json
import re
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "af3_raw"
CASP_DIR = ROOT / "data" / "CASP15"
MSA_DIR = ROOT / "msa"

ZIP_RE = re.compile(r"fold_(t\d+)_(custom|default)_msa\.zip$", re.IGNORECASE)


def read_seq(casp_id: str) -> str:
    fasta = CASP_DIR / f"{casp_id}.fasta"
    return "".join(l.strip() for l in fasta.read_text().splitlines() if not l.startswith(">"))


def import_zip(zip_path: Path):
    m = ZIP_RE.search(zip_path.name)
    if not m:
        print(f"skipping {zip_path.name}: not a fold_<target>_<condition>_msa.zip")
        return
    casp_id, condition = m.group(1).upper(), m.group(2).lower()
    job_name = zip_path.stem  # e.g. fold_t1183_custom_msa
    cond_dir = ROOT / f"af3_{condition}" / casp_id
    out_dir = cond_dir / "server_output"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out_dir)
        # the server writes zip entry times in UTC
        finished_at = max(
            "%04d-%02d-%02dT%02d:%02dZ" % zi.date_time[:5] for zi in zf.infolist()
        )

    request = json.loads((out_dir / f"{job_name}_job_request.json").read_text())[0]
    chain = request["sequences"][0]["proteinChain"]
    if chain["sequence"] != read_seq(casp_id):
        raise SystemExit(f"[{casp_id}/{condition}] job sequence != data/CASP15/{casp_id}.fasta")

    used_msas = sorted((out_dir / "msas").glob("*.a3m"))
    if condition == "custom":
        ours = (MSA_DIR / casp_id / f"{casp_id}_custom.a3m").read_bytes()
        if not any(p.read_bytes() == ours for p in used_msas):
            raise SystemExit(f"[{casp_id}/custom] server MSA differs from msa/{casp_id}/{casp_id}_custom.a3m")

    samples = []
    for summary in sorted(out_dir.glob(f"{job_name}_summary_confidences_*.json")):
        idx = int(summary.stem.rsplit("_", 1)[1])
        samples.append((json.loads(summary.read_text())["ranking_score"], idx))
    best_score, best = max(samples, key=lambda s: (s[0], -s[1]))

    shutil.copy(out_dir / f"{job_name}_model_{best}.cif", cond_dir / f"{casp_id}_model.cif")
    shutil.copy(out_dir / f"{job_name}_full_data_{best}.json", cond_dir / f"{casp_id}_confidence.json")
    shutil.copy(
        out_dir / f"{job_name}_summary_confidences_{best}.json",
        cond_dir / f"{casp_id}_summary_confidences.json",
    )

    template_ids = sorted(
        {p.read_text().splitlines()[0].removeprefix("data_") for p in (out_dir / "templates").glob("*.cif")}
    )
    record = {
        "casp_id": casp_id,
        "condition": condition,
        "job_id": request["name"],
        "submitted_at": None,  # not recorded by the server download
        "finished_at": finished_at,
        "model_version": "AlphaFold Server (AF3), dialect "
        f"{request['dialect']} v{request['version']}",
        "msa_input": (
            str((MSA_DIR / casp_id / f"{casp_id}_custom.a3m").relative_to(ROOT))
            if condition == "custom"
            else "AF3 server default MSA"
        ),
        "model_seeds": request["modelSeeds"],
        "use_structure_template": chain.get("useStructureTemplate"),
        "template_pdb_ids": template_ids,
        "n_samples": len(samples),
        "top_ranked_sample": best,
        "top_ranking_score": best_score,
        "server_msas": [str(p.relative_to(cond_dir)) for p in used_msas],
        "notes": f"imported from af3_raw/{zip_path.name}",
    }
    (cond_dir / "job_metadata.json").write_text(json.dumps(record, indent=2) + "\n")
    print(
        f"[{casp_id}/{condition}] seed={request['modelSeeds'][0]} top sample={best} "
        f"(ranking_score={best_score}) templates={','.join(template_ids)}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip-dir", type=Path, default=RAW_DIR)
    args = parser.parse_args()

    zips = sorted(args.zip_dir.glob("fold_*_msa.zip"))
    if not zips:
        raise SystemExit(f"no fold_*_msa.zip files in {args.zip_dir}")
    for z in zips:
        import_zip(z)


if __name__ == "__main__":
    main()
