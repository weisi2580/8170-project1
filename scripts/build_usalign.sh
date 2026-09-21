#!/usr/bin/env bash
# Build the US-align binary (used for TM-score/RMSD in scripts/evaluate_structures.py)
# from the vendored source in scripts/vendor/ (originally from
# https://github.com/pylelab/USalign, MIT-style academic license -- see
# scripts/vendor/USalign_LICENSE). Vendored so evaluation is reproducible
# without a network dependency; output goes to bin/USalign.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p bin
g++ -O3 -ffast-math -o bin/USalign scripts/vendor/USalign.cpp
echo "built bin/USalign"
bin/USalign -h | head -3
