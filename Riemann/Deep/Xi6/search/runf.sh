#!/bin/bash
# runf.sh T mode spec... : one process per spec, hard alarm T+30 s
T=$1; M=$2; shift 2
for sp in "$@"; do
  out=$(OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 perl -e "alarm $T+30; exec @ARGV" venv/bin/python fam.py "$sp" $M $T 2>&1 | tail -2)
  last=$(echo "$out" | tail -1)
  case "$last" in BUILT*) echo "$last -> TIMEOUT"; echo "{\"spec\": \"$sp\", \"res\": \"TIMEOUT\"}" >> results.jsonl;; *) echo "$last";; esac
done
