# Project 1 Execution Guide for Claude Code

**Course**: CMP_SC 8170 — Project 1
**Topic**: Custom MSA vs. AlphaFold3 default alignment — structure accuracy comparison
**Targets**: T1183 (PDB 8IFX chain B), T1112 (PDB 8ORK), T1122 (PDB 8BBT chain A)

---

## 0. How to use this document

Give this file to Claude Code as the project brief. Each step lists: **Goal**, **Actions**, **Outputs**, and **Manual step?** (whether a human must do part of it, e.g. web-only tools).

**Current folder layout (as of this revision):**
```
project1/
  AI-LOG-template.md
  data/
    CASP15/   T1112.fasta  T1122.fasta  T1183.fasta
    RCSB/     8BBT.cif     8IFX.cif     8ORK.cif
```

---

## Step 1 — Fetch target sequences & reference structures ✅ DONE

**Goal**: For each of the 3 targets, obtain the CASP15 FASTA and the experimental PDB reference.

**Status**: **Completed, including the follow-ups.** FASTA files for T1112, T1122, T1183 are saved in `data/CASP15/` (originally saved as macOS `.textClipping` files by mistake — re-extracted to real `.fasta` text), and reference structures 8BBT, 8IFX, 8ORK are saved (as `.cif`) in `data/RCSB/`.

**Follow-up status**:
- ✅ Chain extraction done: `data/RCSB_chain/8IFX_B.cif` and `data/RCSB_chain/8BBT_A.cif` (8ORK is single-chain, used as-is).
- ✅ `data/metadata.csv` built. Note: `casp_difficulty_class` column is `TBD` — fill in from the official CASP15 classification before using it in the report.
- ✅ `scripts/fetch_targets.py` written (chain extraction + metadata build; `--download` re-fetches CIFs from RCSB).

**Outputs**: `data/CASP15/*.fasta` ✅, `data/RCSB/*.cif` ✅, `data/RCSB_chain/*.cif` ✅, `data/metadata.csv` ✅.

**Manual step?**: Already done manually via the CASP15 and RCSB PDB websites (browser download — no batch API needed for 3 files).

---

## Step 2 — Environment & repo setup ✅ DONE

**Goal**: Reproducible working environment.

**Actions**:
- Extend the existing folder structure with the remaining working directories:
  ```
  project1/
    AI-LOG-template.md
    data/
      CASP15/   T1112.fasta  T1122.fasta  T1183.fasta
      RCSB/     8BBT.cif     8IFX.cif     8ORK.cif
      metadata.csv
    scripts/
    msa/{T1183,T1112,T1122}/
    af3_custom/{T1183,T1112,T1122}/
    af3_default/{T1183,T1112,T1122}/
    eval/{T1183,T1112,T1122}/
    figures/{T1183,T1112,T1122}/
    env.yml
    README.md
  ```
- Create a conda env with: `mmseqs2`, `biopython`, `pandas`, `matplotlib`, `seaborn`, `logomaker` (or `weblogo`), `pymol-open-source`, `tmtools` or `TMalign` binary, `lddt` (OpenStructure) or `US-align`.
- Write `env.yml` and a `Makefile` or shell scripts to reproduce each stage.

**Outputs**: extended `project1/` skeleton, `env.yml`, README with run instructions.

**Manual step?**: No.

**Status**: **Completed.** Full folder skeleton created. No conda available on this machine, so a Python venv (`requirements.txt`) was used instead and documented alongside `env.yml` in `README.md`. TMalign/US-align wasn't available via Homebrew either, so US-align was vendored from source (`scripts/vendor/`) and compiled locally (`scripts/build_usalign.sh` → `bin/USalign`) instead of pulling it from conda.

---

## Step 3 — Build custom MSA ✅ DONE (Step 5 pending for AF3 input)

**Goal**: One high-quality custom MSA per target via MMseqs2 (ColabFold-style search).

**Actions**:
- Use `colabfold_search` or raw `mmseqs2` against **UniRef30** + an environmental database (BFD-lite or MGnify/ColabFold-envdb), reading input from `data/CASP15/<ID>.fasta`.
  ```bash
  colabfold_search data/CASP15/T1183.fasta \
      uniref30_db envdb_db msa/T1183/ \
      --db1 uniref30_2023 --db2 colabfold_envdb_202108
  ```
- Convert result to `.a3m`.
- Write `scripts/filter_msa.py` to filter by: query coverage threshold (e.g. ≥50%), max redundancy (e.g. `hhfilter -id 90`), min identity threshold — parameterize so you can justify choices in the report.
- Validate: query sequence at row 0 matches FASTA exactly, correct A3M formatting, no length mismatches.

**Outputs**: `msa/<ID>/<ID>_custom.a3m` (validated).

**Manual step?**: No, fully scriptable if databases are locally available or accessible via ColabFold's remote MSA server.

**Status**: **Completed for all 3 targets**, using the ColabFold remote MSA server path (`scripts/build_msa.py`, no local databases). `scripts/filter_msa.py` filters + validates (coverage ≥ 0.5, redundancy < 0.9 identity, min identity ≥ 0.2 — all CLI-tunable). Result: T1183 depth 581, T1112 depth 855, **T1122 depth only 2** (its nearest UniRef/env hit is essentially its only homolog — expected for a nudivirus protein, and a useful data point for the hypothesis).

---

## Step 4 — MSA quality control & visualization ✅ DONE

**Goal**: Quantify and visualize MSA depth, coverage, diversity, and conservation.

**Actions**: write `scripts/msa_qc.py` computing, per target:
1. **Depth & Neff/L** — effective sequence count normalized by length.
2. **Query coverage** — per-hit fraction of target residues aligned → coverage distribution plot.
3. **Gap fraction** — per-column gap % → line/heatmap plot.
4. **Identity distribution** — pairwise/query identity histogram.
5. **Conservation/entropy** — per-column Shannon entropy → sequence logo (logomaker/weblogo) + optionally export alignment for Jalview.

**Outputs**: `figures/<ID>/msa_qc_*.png`, a summary CSV comparing all 3 targets (this feeds the "Target comparison table").

**Manual step?**: No.

**Status**: **Completed for all 3 targets** (`scripts/msa_qc.py`). Neff/L: T1183 = 2.63, T1112 = 1.21, T1122 = 0.004. Summary table at `eval/msa_qc_summary.csv`; figures (coverage, gap fraction, identity, entropy, sequence logo) at `figures/<ID>/msa_qc_*.png`.

---

## Step 5 — AlphaFold3 predictions (two conditions)

**Goal**: Run each target twice — once with the custom A3M, once with AF3's default alignment — keeping every other setting identical.

**Actions**:
- **Condition A (custom)**: Upload `msa/<ID>/<ID>_custom.a3m` + target sequence to the AlphaFold3 server; record job ID, submission timestamp, and settings in `af3_custom/<ID>/job_metadata.json`.
- **Condition B (default)**: Submit the same raw sequence, let AF3 build its own MSA; record the same metadata fields in `af3_default/<ID>/job_metadata.json`.
- Download for both: top-ranked model (`.cif`/`.pdb`), confidence JSON (pLDDT, PAE).
- Write `scripts/log_af3_job.py` to standardize the metadata JSON schema across all 6 runs (3 targets × 2 conditions).

**Outputs**: `af3_custom/<ID>/<ID>_model.cif`, `af3_custom/<ID>/<ID>_confidence.json`, same for `af3_default/`, plus `job_metadata.json` per run.

**Manual step?**: **Yes** — the public AlphaFold3 server is a browser-based tool with no general public batch API for this use case. Claude Code can prepare the exact input files and a checklist, but a human must perform the actual upload/submission/download for each of the 6 jobs, then hand the downloaded files back to Claude Code for the next steps.

**Status**: **Inputs prepared, submissions pending.** `scripts/log_af3_job.py --prep` wrote `input_sequence.fasta` + `SUBMIT_CHECKLIST.txt` into all 6 of `af3_custom/<ID>/` and `af3_default/<ID>/`. Still needed: the human submits all 6 jobs at alphafoldserver.com, logs each with `--log`, and drops the downloaded model + confidence JSON into the matching folder.

---

## Step 6 — Structure evaluation

**Goal**: Score both conditions against the experimental PDB reference.

**Actions**: write `scripts/evaluate_structures.py`:
1. **Global accuracy**: TM-score and RMSD via `TMalign`/`US-align`, restricted to the CASP "evaluation unit" residue range per target, using the chain-filtered reference structures from `data/RCSB/` (see Step 1 follow-up).
2. **Local accuracy**: per-residue lDDT (or CA-RMSD per residue) comparing prediction vs. reference.
3. **Model confidence**: parse pLDDT and PAE directly from AF3's confidence JSON (no reference needed).
4. Aggregate into one row per (target, condition): TM-score, RMSD, mean pLDDT, mean lDDT.

**Outputs**: `eval/scores_summary.csv` (6 rows: 3 targets × 2 conditions).

**Manual step?**: No, once model files are downloaded.

**Status**: **Script written and dry-run tested** (`scripts/evaluate_structures.py`, using vendored/compiled US-align instead of a separate TMalign+lddt toolchain — see Step 2). Currently produces no rows since no AF3 models exist yet; will populate `eval/scores_summary.csv` automatically once Step 5 lands models in `af3_custom/*/` and `af3_default/*/`.

---

## Step 7 — Visualization & final comparison

**Goal**: Produce all figures and tables needed for the report/slides.

**Actions**:
- PyMOL script (`scripts/overlay.pml` or PyMOL Python API) to overlay: reference vs. custom-MSA model vs. default-MSA model, colored by pLDDT or by RMSD deviation — export as image per target.
- Per-residue confidence plot (pLDDT along sequence) for both conditions, overlaid.
- Metric comparison bar/scatter plot: TM-score and RMSD, custom vs. default, per target.
- Final "target summary table": all metrics + CASP difficulty class + expected homolog depth, to directly support/refute the working hypothesis.
- Write a short `RESULTS.md` stating, per target, whether the custom MSA improved accuracy, and an overall conclusion across the 3 targets.

**Outputs**: `figures/<ID>/overlay.png`, `figures/<ID>/plddt.png`, `figures/metric_comparison.png`, `eval/target_summary_table.csv`, `RESULTS.md`.

**Manual step?**: No (overlay.pml needs PyMOL installed, which isn't set up on this machine yet — see env.yml).

**Status**: **Scripts written and dry-run tested** (`scripts/plot_comparison.py`, `scripts/overlay.pml`, `scripts/write_results.py`). `eval/target_summary_table.csv` and a placeholder `RESULTS.md` already generate correctly from the Step 1–4 data with AF3 columns blank; they'll fill in once Step 5/6 provide real scores.

---

## Step 8 — AI-usage log (for the "How we used AI" slide)

**Goal**: Keep a running log of exactly what Claude Code did at each step, for transparency in the report.

**Actions**: Use/extend the existing `AI-LOG-template.md`. Append one line per script/action: date, step, what was generated/run, what a human reviewed or changed.

---

## Suggested run order

1 ✅, 2 ✅, 3 ✅, 4 ✅ — all completed in one Claude Code session, with real data (see per-step Status notes above).
5 requires a manual pause for the 6 AF3 web submissions — **currently blocking** (inputs/checklists are ready in `af3_custom/*/` and `af3_default/*/`).
6→7→8 are scripted and dry-run tested against empty AF3 output; they'll resume automatically and produce real numbers once the model files are back on disk.