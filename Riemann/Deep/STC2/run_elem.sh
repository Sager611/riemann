#!/bin/bash
# sup F for all six kernels, modes N and S, [1,1e4] (the paper's T0) by intervals, Lemma 4.7 beyond; one process at a time under memcap
cd /Users/adriansager/playground/riemann; export OMP_NUM_THREADS=1
for kv in L5:0.20176 L20:0.10823 L60:0.09473 L150:0.09349 C400:0.08399 D6000:0.08382; do
  for mode in N S; do
    bash .scratch/memcap.sh 290 900 uv run --no-project --with python-flint python3 Riemann/Deep/STC2/elem_arb.py ${kv%%:*} ${kv##*:} $mode 1e4
  done
done
