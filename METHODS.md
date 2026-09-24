# Methods: every setting, and why the custom MSA lost

This page collects every setting of the custom MSA and the AF3 runs, lists every way the two conditions differ, and records the evidence for why the custom-MSA models were less accurate. Numbers come from `msa/`, `af3_*/*/job_metadata.json` and `eval/*.csv`; the template analysis is `scripts/analyze_templates.py` → `eval/template_analysis.csv`, `eval/template_summary.csv`.

## 1. Custom MSA: what we set

| Stage | Setting | Where |
|---|---|---|
| Search service | Public ColabFold MMseqs2 API, `https://api.colabfold.com/ticket/msa` | `scripts/build_msa.py` |
| Search mode | `mode=env`: UniRef search (ColabFold's UniRef30 profile search, hits expanded to UniRef100 members) **+** ColabFold environmental DB (BFD, MGnify, MetaEuk, SMAG). No paired MSA requested. | `build_msa.py` |
| Search parameters | ColabFold server defaults, not set by us. From the server's own `msa/<ID>/raw/msa.sh`: 3 iterations, `-e 0.1`, `--max-seqs 10000`, realign `-e 10`, server-side diversity filter `--diff 3000 --qid 0,0.2,…,1.0 --max-seq-id 0.95` | `msa/<ID>/raw/msa.sh` |
| Templates from ColabFold | Downloaded (`pdb70.m8`) but **not used**. AF3 ran its own template search in both conditions. | — |
| Merge | `uniref.a3m` + `bfd.mgnify30.metaeuk30.smag30.a3m`, query row first, exact-duplicate aligned rows removed → `msa/<ID>/<ID>_raw.a3m` | `build_msa.py` |
| Filter 1: query coverage | keep a hit only if it aligns to **≥ 50%** of query positions (`--min-coverage 0.5`) | `scripts/filter_msa.py` |
| Filter 2: identity to query | keep only if **≥ 20%** identical to the query over columns where both have a residue (`--min-identity 0.2`) | `filter_msa.py` |
| Filter 3: redundancy | drop a hit **≥ 90%** identical to any hit already kept, greedy in file order: UniRef hits first, then environmental (`--max-identity 0.9`) | `filter_msa.py` |
| Validation | row 0 == CASP15 FASTA exactly; equal match-column length for every row; A3M alphabet only | `filter_msa.py` |
| Upload | `msa/<ID>/<ID>_custom.a3m` uploaded to the AF3 Server as chain A's unpaired MSA. On the server it **replaces** AF3's own MSA; there is no paired MSA. Byte-identity with the uploaded file is checked by `scripts/import_af3_results.py`. | manual + import check |

### How many sequences survived each step

| Target | UniRef hits | Env hits | Merged (raw) | Removed: coverage < 50% | Removed: identity < 20% | Removed: ≥ 90% redundant | **Final custom MSA** | AF3 default MSA (unpaired + paired) |
|---|---|---|---|---|---|---|---|---|
| T1183 | 842 | 1,630 | 2,441 | **1,813** | 0 | 48 | **581** | 8,032 |
| T1112 | 399 | 2,149 | 2,546 | **1,473** | 0 | 219 | **855** | 1,967 |
| T1122 | 1 | 0 | 1 | 0 | 0 | 0 | **2** | 3 |

Hit counts exclude the query row; the two MSA-size columns include it. Before filtering, the T1112 search had **more** sequences than AF3's MSA (2,547 vs 1,967). The coverage filter removed 74% (T1183) and 58% (T1112) of all hits. Half of the raw hits cover less than ~45% of the query.

## 2. AF3 Server runs: what was the same

| Setting | Both conditions |
|---|---|
| Service | AlphaFold Server (AF3), job dialect `alphafoldserver` v3, run 2026-09-23 |
| Input | 1 protein chain, sequence == CASP15 FASTA (checked), 1 copy, no ligands |
| Templates | **on** (`useStructureTemplate: true`); AF3 searched its own templates in every job, up to 4 used |
| Samples | 5 per job; the top-ranked (highest `ranking_score`) is the reported model |
| Seeds | one per job, chosen by the server (not settable in the UI we used): custom 311040490 / 1145570114 / 517604978; default 587090016 / 2099461411 / 775453839 (T1183 / T1112 / T1122) |

## 3. Every way the two conditions differ

| | Custom | AF3 default |
|---|---|---|
| Search tool / databases | MMseqs2 (ColabFold): UniRef + BFD/MGnify/MetaEuk/SMAG | AF3's pipeline (jackhmmer-type search; hits in the returned MSA are UniRef90, MGnify and UniProt entries) |
| Filtering | our coverage ≥ 50%, identity ≥ 20%, redundancy < 90% | AF3's own; keeps partial-coverage fragments |
| Paired MSA | none | yes (UniProt; 2,999 / 507 / 2 rows) — for a monomer these rows are still used as extra MSA rows |
| Template hits and their alignment | found by AF3 from our MSA | found by AF3 from its own MSA |
| Random seed | different per job | different per job |

## 4. Why was the custom model worse? The evidence

### Hypothesis A — "the default run had templates and the custom run did not": **ruled out**

Both conditions ran with templates on and both used 4 template hits (`job_metadata.json → template_pdb_ids`). What differs is *which* templates and *how they are aligned*: AF3 derives the template search from the input MSA. `scripts/analyze_templates.py` scores every template as if it were a partial model of the target, against the experimental structure:

| Target | Condition | Templates | Coverage of target (union) | Best single template, TM vs experiment |
|---|---|---|---|---|
| T1183 | custom | 3zet 4wq4 3r6m 1okj | 99.5% | 0.641 |
| T1183 | default | 3r6m 3zet 4wq4 1okj | 99.5% | 0.691 |
| T1112 | custom | 2hf8 2wsm 1vl8 1vl8 | 34% | 0.161 |
| T1112 | default | 2hf8 2wsm 4lps 3dm5 | 35% | 0.163 |
| T1122 | custom | 2pyb 3f9i 3f9i 4ayd | 77% | 0.114 |
| T1122 | default | 2pyb 4oyu 3f9i 3f9i | 80% | 0.114 |

- **T1183:** same four PDB entries, but aligned differently. The default alignments put the templates closer to the experiment (best TM 0.69 vs 0.64), yet both final models reach TM ≈ 0.98, so this difference did not matter.
- **T1112:** templates cover only ~35% of the chain in either condition and are equally weak (TM ≤ 0.16). The custom-only template 1vl8 sits exactly on the region where the custom model is worse (residues 342–411). But the models did **not** copy it: over those residues the template is 4.4–6.9 Å from the experiment, while *both* models are 1.3–1.8 Å from the experiment and equally far (~4.7–6.6 Å) from the template.
- **T1122:** every template is poor (TM ≤ 0.11, 11–17 Å RMSD over its aligned part) and neither model follows them (models are further from the templates than from the experiment).

So templates do not explain the gap. They add almost no correct information on T1112 and T1122, and the models are not following them.

### Hypothesis B — "the custom MSA lost inter-domain signal": **supported for T1112**

The T1112 error is not a local misfold. Superposing each model on the experiment by region:

| T1112 | Core (1–339), own fit | C-terminal region (340–460), own fit | C-terminal after superposing on the core |
|---|---|---|---|
| custom | 1.90 Å | 2.28 Å | **6.17 Å** |
| default | 1.60 Å | 2.02 Å | **4.27 Å** |

Both models fold the C-terminal region nearly as well locally. The custom model **places it at the wrong angle** relative to the core. Domain orientation comes from co-evolution between the two domains, which needs many full-length homologs. That is exactly what the ≥ 50% coverage filter thinned out (1,473 hits removed) and what AF3's paired rows add. Templates cannot supply it (none spans both domains).

### Hypothesis C — "random seed / sampling": **cannot be excluded for T1122**

On T1122 both MSAs are effectively single-sequence and the templates are uninformative, so there is no systematic input difference left that could explain a 0.09 TM gap. The 5 samples within each job agree closely (sd ≤ 0.01), but samples share the job's seed, so seed-to-seed variation is not measured. The honest reading for T1122: the difference is most likely run-to-run variation, not a property of the MSA.

## 5. Experiments that would settle it (need AF3 Server submissions)

1. **Unfiltered custom MSA:** upload `msa/<ID>/<ID>_raw.a3m` (no coverage filter) for T1112. If TM approaches 0.93, the coverage filter is the cause.
2. **Templates off, both conditions:** removes template differences entirely.
3. **Repeat each job with 2–3 more seeds:** measures seed-to-seed variance, which decides T1122.

## 6. Evaluation settings (for reference)

- Scoring region: CASP15 evaluation units (T1183 1–195, T1112 1–460, T1122 4–237).
- TM-score and Cα RMSD: US-align `-TMscore 1` (residue-index correspondence) on Cα files written by `scripts/evaluate_structures.py` (MSE → MET), TM normalized by the reference length.
- Cα-lDDT: inclusion radius 15 Å, thresholds 0.5/1/2/4 Å, mean over residues.
- MSA Neff: sequences weighted by 1/(number within 80% identity), divided by L.
