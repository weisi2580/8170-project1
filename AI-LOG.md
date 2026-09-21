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

## Example entry (delete before first commit)

### [2026-10-14] — First pass at contributor-activity view
- **Tool/model:** Claude Sonnet 4.6
- **What I asked for:** a multi-series time chart of monthly commit counts for the top 8 contributors, from my cleaned dataframe
- **What it produced:** working Plotly code using a stacked area chart with a qualitative palette
- **Accepted:** the data-wrangling (groupby/resample) after checking totals against a manual count for two months
- **Rejected/modified:** the stacked area encoding. My task is comparing individual contributors' trends; stacking puts all but the bottom series on a shifting baseline, which destroys position-along-common-scale judgment (Cleveland). Switched to small multiples with shared y-axis.
- **Verification:** re-ran against raw data; spot-checked contributor #3's spike against the actual commit log for that month.

---

### [2026-09-21] — Repo setup + Steps 1–4 (chain extraction, MSA build/filter/QC) for all 3 targets
- **Tool/model:** Claude Sonnet 5 (Claude Code CLI)
- **What I asked for:** execute project_plan.md end to end — build the repo structure, README, and run the pipeline as far as it can go without the manual AF3 web step.
- **What it produced:** discovered the 3 "FASTA" files were actually broken macOS `.textClipping` files and recovered the real sequences from them; wrote `scripts/fetch_targets.py` (chain extraction for 8IFX→B and 8BBT→A via Biopython, `data/metadata.csv`); wrote `scripts/build_msa.py` implementing the ColabFold remote MMseqs2 search API (no local databases) and ran it live for T1183/T1112/T1122; wrote `scripts/filter_msa.py` (coverage/redundancy/identity filtering + validation) and `scripts/msa_qc.py` (Neff/L, coverage/identity/entropy plots, sequence logo) and ran both on real data; since Homebrew had no TMalign/US-align formula, vendored and compiled US-align from source (`scripts/vendor/`, `scripts/build_usalign.sh`) instead; wrote `scripts/evaluate_structures.py`, `scripts/plot_comparison.py`, `scripts/write_results.py`, `scripts/overlay.pml`, `scripts/log_af3_job.py` (AF3 input prep for the manual step) and dry-ran all of them against the current (AF3-less) state; wrote `README.md`, `env.yml`, `requirements.txt`, `Makefile`, `.gitignore`.
- **Accepted:** all of the above as-is after checking outputs — e.g. confirmed `data/RCSB_chain/8IFX_B.cif` actually dropped from 4479 to 1738 atom records, confirmed `msa/T1122/T1122_custom.a3m` genuinely only has 2 rows (not a bug — verified the raw server response), confirmed `scripts/fetch_targets.py`'s metadata.csv matches the plan's target→PDB→chain mapping exactly.
- **Rejected/modified:** first version of the identity metric used two different denominators in `filter_msa.py` (hit-only non-gap count) vs `msa_qc.py` (union of non-gap positions in either sequence) — caught because a summary row showed mean_identity (0.18) below the 0.2 filter floor, which is impossible if both used the same definition. Fixed both to use the same "both sequences non-gap" convention (statistical/engineering consistency) and re-ran Steps 3b/4 on all 3 targets.
- **Verification:** re-ran the full filter+QC pipeline after the identity-metric fix and confirmed mean_identity ≥ min_identity threshold for all targets; spot-checked raw MSA depths against `wc -l` on the actual `.a3m` files; ran `scripts/evaluate_structures.py`/`plot_comparison.py`/`write_results.py` against the current no-AF3-data state to confirm they fail gracefully (print "no model found yet" / write a placeholder RESULTS.md) rather than crashing, before handing off the manual AF3 step to the user.
