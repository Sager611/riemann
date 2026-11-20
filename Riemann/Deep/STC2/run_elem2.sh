#!/bin/bash
# cross-check: sup F on [1,1e8] by intervals (second-order forms beyond 1e4), Lemma 4.7 only for T >= 1e8; one process at a time
cd /Users/adriansager/playground/riemann; export OMP_NUM_THREADS=1
for kv in L5:0.20176 L20:0.10823 L60:0.09473 L150:0.09349 C400:0.08399 D6000:0.08382; do
  for mode in N S; do
    bash .scratch/memcap.sh 290 900 uv run --no-project --with python-flint python3 Riemann/Deep/STC2/elem_arb2.py ${kv%%:*} ${kv##*:} $mode 1e8
  done
done
