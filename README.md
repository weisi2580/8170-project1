# Custom MSA vs. AlphaFold3 Default Alignment

**Course:** CMP_SC 8170 — Project 1
**Question:** does a custom MMseqs2/ColabFold multiple sequence alignment (MSA) improve AlphaFold3 structure prediction accuracy over AF3's own default alignment?
**Targets:** three CASP15 targets spanning a range of expected homolog depth:

| CASP ID | PDB (chain) | Description |
|---|---|---|
| T1183 | 8IFX (chain B) | *Aquifex aeolicus* TsaB, complexed with TsaD |
| T1112 | 8ORK | *Methanothermus fervidus* protein |
| T1122 | 8BBT (chain A) | *Tipula oleracea* nudivirus protein (expected shallow MSA) |

For each target, AF3 is run twice — once seeded with a custom MSA built via the ColabFold remote MMseqs2 search server, once letting AF3 build its own default alignment — with every other setting held identical. Both predictions are then scored against the experimental PDB structure with US-align (TM-score, RMSD) and compared against AF3's own pLDDT confidence.

## Repository structure

```
project1/
├── data/
│   ├── CASP15/            target FASTAs (T1112, T1122, T1183)
│   ├── RCSB/               full experimental reference structures (8BBT, 8IFX, 8ORK)
│   ├── RCSB_chain/         single-chain-filtered references used for scoring
│   │                       (8IFX_B.cif, 8BBT_A.cif; 8ORK is already single-chain)
│   └── metadata.csv        per-target CASP ID, length, PDB ID, chain, download date, difficulty class
├── msa/<ID>/
│   ├── raw/                 unmerged per-database hits from the ColabFold API (gitignored, regenerable)
│   ├── <ID>_raw.a3m         merged UniRef100 + environmental hits, before filtering
│   └── <ID>_custom.a3m      final filtered/validated custom MSA fed to AF3
├── af3_custom/<ID>/         AF3 "custom MSA" condition: input, submission checklist, downloaded model + confidence JSON, job_metadata.json
├── af3_default/<ID>/        AF3 "default MSA" condition: same, but AF3 builds its own alignment
├── eval/                    scores_summary.csv, msa_qc_summary.csv, target_summary_table.csv
├── figures/<ID>/             all QC and comparison plots
├── scripts/                  every pipeline stage (see below), plus vendor/ (US-align source)
├── bin/USalign               compiled structure-alignment binary (gitignored, build with scripts/build_usalign.sh)
├── env.yml                   conda environment spec
├── requirements.txt          pip freeze of the same deps, for the venv path actually used to develop this repo
└── Makefile                  `make <stage>` wrapper around the scripts below
```

## Setup

This was developed without conda available, using a plain Python venv; either path works.

**Option A — conda (as specified in the original project plan):**
```bash
conda env create -f env.yml
conda activate casp15-msa-af3
```

**Option B — venv (what this repo was actually built with):**
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Either way, also build the structure-alignment tool (vendored source, no network needed):
```bash
scripts/build_usalign.sh    # writes bin/USalign
```

## How to run

The pipeline mirrors `project_plan.md`'s steps 1–8. Steps 1–4 are fully automatic and have already been run once (outputs are committed); steps 6–8 will re-run automatically once the AF3 models from step 5 are downloaded.

| Stage | Command | Notes |
|---|---|---|
| 1. Fetch/verify targets + chain-filter references | `python scripts/fetch_targets.py` | Add `--download` to re-fetch CIFs from RCSB |
| 2. (this file / env setup) | — | done |
| 3. Build custom MSA | `python scripts/build_msa.py [--target T1183]` | Talks to the public ColabFold MMseqs2 API (`api.colabfold.com`) — no local database needed. Be polite: don't run targets in parallel. |
| 3b. Filter MSA | `python scripts/filter_msa.py --target T1183` | Coverage / redundancy / min-identity thresholds are CLI flags; validates row 0 == query FASTA, uniform alignment length, legal A3M alphabet |
| 4. MSA QC + figures | `python scripts/msa_qc.py --target T1183` | Depth, Neff/L, coverage & identity distributions, per-column gap fraction & entropy, sequence logo |
| 5. AF3 predictions (**manual**) | `python scripts/log_af3_job.py --prep` then submit by hand | See "Step 5 is manual" below |
| 6. Evaluate structures | `python scripts/evaluate_structures.py` | TM-score/RMSD via US-align against `data/RCSB_chain/`; mean pLDDT from AF3's confidence JSON |
| 7. Comparison figures + table | `python scripts/plot_comparison.py` | `figures/metric_comparison.png`, per-target `plddt.png`, `eval/target_summary_table.csv` |
| 7b. Overlay figures | `pymol -cq scripts/overlay.pml -- T1183` | Requires PyMOL (`pymol-open-source` in env.yml); one call per target |
| 7c. Results writeup | `python scripts/write_results.py` | Generates `RESULTS.md` from the summary table |

Or via `make`:
```bash
make setup usalign
make all           # stages 1, 3, 3b, 4, and AF3 input prep
# ... manual AF3 step (see below) ...
make evaluate compare results
```

### Step 5 is manual

The public AlphaFold3 server ([alphafoldserver.com](https://alphafoldserver.com)) has no batch API for this use case. `scripts/log_af3_job.py --prep` writes the exact sequence and a checklist into each `af3_custom/<ID>/SUBMIT_CHECKLIST.txt` and `af3_default/<ID>/SUBMIT_CHECKLIST.txt` — follow those for all 6 jobs (3 targets × 2 conditions), then:
1. record each job with `python scripts/log_af3_job.py --log --target <ID> --condition <custom|default> --job-id <ID> --submitted-at <ISO8601> --model-version <version>`
2. drop the downloaded top-ranked model (`.cif`/`.pdb`) and confidence JSON into the same folder.

## Current status

- **Steps 1–4: done**, with real data — custom MSAs for all 3 targets were built live against the ColabFold remote server (`msa/<ID>/<ID>_custom.a3m`), filtered and QC'd (`eval/msa_qc_summary.csv`, `figures/<ID>/msa_qc_*.png`). Reference structures are chain-filtered (`data/RCSB_chain/`) and `data/metadata.csv` is built.
  - As expected, T1122 (a nudivirus protein) got a very shallow MSA (depth 2, Neff/L ≈ 0.004) versus T1183 and T1112 (depth in the hundreds, Neff/L ≈ 1.2–2.6) — directly relevant to the working hypothesis that custom-MSA gains should be largest for shallow/hard targets.
- **`data/metadata.csv`'s `casp_difficulty_class` column is `TBD`** — fill it in from the official CASP15 classification at predictioncenter.org before using the table in the report; this repo does not guess it.
- **Step 5 (AF3 predictions) has not been run yet** — inputs and checklists are prepared in `af3_custom/*/` and `af3_default/*/`, waiting on the 6 manual browser submissions.
- **Steps 6–8 are wired up and tested against empty inputs** (they report "no model found yet" cleanly); they'll produce real numbers as soon as the AF3 outputs land. See `RESULTS.md` for the current (placeholder) state.

## Notes on tooling substitutions

- **MSA search**: rather than a local `mmseqs2` install against multi-hundred-GB UniRef30/BFD databases, `scripts/build_msa.py` uses the same public ColabFold MMseqs2 search API (`api.colabfold.com`) that the ColabFold notebook uses — this is the "remote MSA server" path the project plan explicitly allows for Step 3.
- **Structure alignment**: TMalign/US-align isn't in Homebrew and wasn't installed locally, so `scripts/vendor/` carries the US-align source ([github.com/pylelab/USalign](https://github.com/pylelab/USalign)) and `scripts/build_usalign.sh` compiles it with a plain C++ compiler — no conda/database dependency.

## AI usage log

See `AI-LOG.md` for a running log of AI-assisted work (kept outside the public repo — gitignored per project convention; not part of this GitHub repository).
