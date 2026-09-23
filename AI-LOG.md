# AI-LOG

**Course:** CS 8630 | **Student:** _name_ | **Project/Assignment:** _which_

> Commit an entry in the same commit as the work it describes. Entries are graded on
> judgment quality, not volume. Summaries, not transcripts. One entry per meaningful
> AI interaction or work session — not per prompt.

---

## Entry template

### [YYYY-MM-DD] — _short description of the work session_
- **Tool/model:** (e.g., Claude Sonnet 4.6, GPT-5, Copilot in VS Code)
- **What I asked for:** one or two sentences summarizing intent, not the verbatim prompt
- **What it produced:** one or two sentences
- **Accepted:** what you kept, and why it was right
- **Rejected/modified:** what you threw out or changed, and the *principle* that justified it
  (perceptual, grammatical, statistical, or engineering grounds — name it)
- **Verification:** how you checked correctness (ran it, compared to docs, tested edge case,
  cross-checked the statistic)

---

### [2026-09-21] — Repo setup + Steps 1–4 (chain extraction, MSA build/filter/QC) for all 3 targets
- **Tool/model:** Claude Sonnet 5 (Claude Code CLI)
- **What I asked for:** execute project_plan.md end to end — build the repo structure, README, and run the pipeline as far as it can go without the manual AF3 web step.
- **What it produced:** discovered the 3 "FASTA" files were actually broken macOS `.textClipping` files and recovered the real sequences from them; wrote `scripts/fetch_targets.py` (chain extraction for 8IFX→B and 8BBT→A via Biopython, `data/metadata.csv`); wrote `scripts/build_msa.py` implementing the ColabFold remote MMseqs2 search API (no local databases) and ran it live for T1183/T1112/T1122; wrote `scripts/filter_msa.py` (coverage/redundancy/identity filtering + validation) and `scripts/msa_qc.py` (Neff/L, coverage/identity/entropy plots, sequence logo) and ran both on real data; since Homebrew had no TMalign/US-align formula, vendored and compiled US-align from source (`scripts/vendor/`, `scripts/build_usalign.sh`) instead; wrote `scripts/evaluate_structures.py`, `scripts/plot_comparison.py`, `scripts/write_results.py`, `scripts/overlay.pml`, `scripts/log_af3_job.py` (AF3 input prep for the manual step) and dry-ran all of them against the current (AF3-less) state; wrote `README.md`, `env.yml`, `requirements.txt`, `Makefile`, `.gitignore`.
- **Accepted:** all of the above as-is after checking outputs — e.g. confirmed `data/RCSB_chain/8IFX_B.cif` actually dropped from 4479 to 1738 atom records, confirmed `msa/T1122/T1122_custom.a3m` genuinely only has 2 rows (not a bug — verified the raw server response), confirmed `scripts/fetch_targets.py`'s metadata.csv matches the plan's target→PDB→chain mapping exactly.
- **Rejected/modified:** first version of the identity metric used two different denominators in `filter_msa.py` (hit-only non-gap count) vs `msa_qc.py` (union of non-gap positions in either sequence) — caught because a summary row showed mean_identity (0.18) below the 0.2 filter floor, which is impossible if both used the same definition. Fixed both to use the same "both sequences non-gap" convention (statistical/engineering consistency) and re-ran Steps 3b/4 on all 3 targets.
- **Verification:** re-ran the full filter+QC pipeline after the identity-metric fix and confirmed mean_identity ≥ min_identity threshold for all targets; spot-checked raw MSA depths against `wc -l` on the actual `.a3m` files; ran `scripts/evaluate_structures.py`/`plot_comparison.py`/`write_results.py` against the current no-AF3-data state to confirm they fail gracefully (print "no model found yet" / write a placeholder RESULTS.md) rather than crashing, before handing off the manual AF3 step to the user.

### [2026-09-23] — Steps 5–7: import AF3 results, evaluate 30 models, figures, RESULTS.md
- **Tool/model:** Claude Opus 5.5 (Claude Code CLI)
- **What I asked for:** check whether my 6 AlphaFold Server downloads (3 targets × custom/default MSA) were what the pipeline needed, then finish the project and work out how to present it.
- **What it produced:** checked the zips (job sequence == FASTA; uploaded MSA byte-identical to `msa/<ID>/<ID>_custom.a3m`; no template is a reference structure, all ≤ 2015). Moved them to `af3_raw/` and wrote `scripts/import_af3_results.py`. Rewrote `scripts/evaluate_structures.py` to score all 5 samples per run on the CASP15 evaluation unit (US-align `-TMscore 1` TM/RMSD, numpy Cα-lDDT, pLDDT/PAE/pTM). Added `scripts/compare_msas.py` (our MSA vs the MSA AF3 built), a new `scripts/plot_comparison.py` (metric comparison, MSA-depth-vs-TM, pLDDT-vs-lDDT calibration, per-residue coverage/pLDDT/error, PAE), `scripts/render_overlays.py` (PyMOL pip wheel; replaces `overlay.pml`), and a data-driven `scripts/write_results.py` → `RESULTS.md`. Filled the CASP15 difficulty class + evaluation units from predictioncenter.org, updated README/Makefile (`make analysis`), and published an interactive results page (3D overlay + per-residue tracks) for the presentation.
- **Accepted:** the headline result, after checking it is robust. The custom MSA never beat AF3's default (ΔTM −0.003 / −0.038 / −0.087 for T1183 / T1112 / T1122; for T1112 and T1122 the 5-sample ranges don't overlap), so the working hypothesis is not supported. Also accepted the explanation that AF3's server MSA is deeper (8032 vs 581, 1967 vs 855 sequences) and that for T1122 both MSAs are effectively single-sequence.
- **Rejected/modified:** (1) The first T1112 scores were wrong. US-align silently drops the 12 selenomethionine (MSE, HETATM) residues in 8ORK (it read 448 of 460), so TM was normalized over the wrong residue set. Now it is fed clean Cα files and asserts the residue count (measurement validity). (2) Removed claims in the draft write-up that went beyond the data. "Default better everywhere" on T1112 became "equal or better lDDT at ~76% of residues". "Template search follows the MSA" became "consistent with". "Unrelated fragment" became "37-residue fragment, 35% identity" (don't state causes that weren't tested). (3) Changed the Neff comparison to count AF3's paired + unpaired MSA, since a monomer sees both (a like-for-like comparison). (4) Fixed chart issues: labels colliding with titles, a log axis on a depth-2 MSA, and clipped/grey PyMOL renders.
- **Verification:** PyMOL's independent least-squares fit reproduces every US-align RMSD exactly (0.70/0.63, 3.00/2.28, 18.47/14.17 Å). Deleted all derived outputs and ran `make analysis` from the zips: `RESULTS.md` and `target_summary_table.csv` came out byte-identical. Cross-checked the CASP classes/evaluation units against the raw predictioncenter.org HTML. Spot-checked the per-region claims (T1112 340–425, T1122 helical segments) against `eval/<ID>/per_residue_*.csv`.
