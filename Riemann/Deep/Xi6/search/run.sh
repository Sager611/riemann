#!/bin/bash
# run.sh T name... : one process per instance, hard alarm T seconds
T=$1; shift
for nm in "$@"; do
  out=$(OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 perl -e "alarm $T; exec @ARGV" venv/bin/python batch.py 100000 "$nm" 2>/dev/null | tail -1)
  if [ -z "$out" ]; then echo "{\"name\": \"$nm\", \"res\": \"TIMEOUT>$T\"}" | tee -a results.jsonl; else echo "$out"; fi
done
