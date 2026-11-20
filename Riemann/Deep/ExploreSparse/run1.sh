#!/bin/bash
# explore61 run1: P = R + Q split, zeta, head theta in (1/K0, 1) on the q-grid Q/k.  Output: logs/s1_scan.log
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
for args in "$@"; do
  /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1500 uv run --no-project --with numpy --with scipy python3 s1_split.py $args 2>&1 \
   | rg 'N=|memcap|maximum resident|peak memory' | tr '\n' ' ' | sed 's/$/\n/' >> logs/s1_scan.log
  tail -1 logs/s1_scan.log | cut -c1-260
done
