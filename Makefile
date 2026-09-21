TARGETS = T1183 T1112 T1122
PY = python3

.PHONY: setup usalign fetch msa filter qc af3-prep evaluate compare results all clean

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
# and af3_default/<ID>/SUBMIT_CHECKLIST.txt. Run the targets below after the
# 6 AF3 model + confidence files have been downloaded into place.

evaluate:
	$(PY) scripts/evaluate_structures.py

compare:
	$(PY) scripts/plot_comparison.py

results:
	$(PY) scripts/write_results.py

# stages 1-4: fully automatic, no AF3 outputs needed yet
all: fetch usalign msa filter qc af3-prep

clean:
	rm -rf msa/*/raw msa/*/*_raw.a3m msa/*/*_custom.a3m
	rm -rf figures/*/*.png
	rm -f eval/*.csv RESULTS.md
