TARGETS = T1183 T1112 T1122
PY = python3

.PHONY: setup usalign fetch msa filter qc af3-prep import evaluate msa-compare compare overlays results analysis all clean

setup:
	$(PY) -m venv .venv
	. .venv/bin/activate && pip install -q --upgrade pip && pip install -q -r requirements.txt

usalign:
	scripts/build_usalign.sh

fetch:
	$(PY) scripts/fetch_targets.py

msa:
	$(PY) scripts/build_msa.py

filter:
	for t in $(TARGETS); do $(PY) scripts/filter_msa.py --target $$t; done

qc:
	for t in $(TARGETS); do $(PY) scripts/msa_qc.py --target $$t; done

af3-prep:
	$(PY) scripts/log_af3_job.py --prep

# Step 5 (AF3 submissions) is manual -- see af3_custom/<ID>/SUBMIT_CHECKLIST.txt
# and af3_default/<ID>/SUBMIT_CHECKLIST.txt. Put the 6 downloaded
# fold_<target>_<condition>_msa.zip files in af3_raw/, then run `make analysis`.

import:
	$(PY) scripts/import_af3_results.py

evaluate:
	$(PY) scripts/evaluate_structures.py

msa-compare:
	$(PY) scripts/compare_msas.py

compare:
	$(PY) scripts/plot_comparison.py

overlays:
	$(PY) scripts/render_overlays.py

results:
	$(PY) scripts/write_results.py

# stages 5b-7: everything downstream of the AF3 zips
analysis: import evaluate msa-compare compare overlays results

# stages 1-4: fully automatic, no AF3 outputs needed yet
all: fetch usalign msa filter qc af3-prep

clean:
	rm -rf msa/*/raw msa/*/*_raw.a3m msa/*/*_custom.a3m
	rm -rf figures/*/*.png
	rm -rf af3_custom/*/server_output af3_default/*/server_output
	rm -f eval/*.csv eval/*/*.csv RESULTS.md
