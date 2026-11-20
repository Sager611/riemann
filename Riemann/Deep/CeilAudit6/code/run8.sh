#!/bin/bash
# ceilaudit6 run8: one window model containing 8 pockets (S = {inf,2,3,5,7}, tc = 23150, tau1 = 1200): count of angles vs pockets.
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy python3"
echo "== depth check (tc 35053.2, tmax 100): depth 6 vs 16"
for dp in 6 16; do $UV wc_copy.py 2,3,5,7 35053.2 100 $dp | cut -c1-330; done
echo "== multi-pocket window"
python3 pk_in.py logs/pk_2357_20000.json 23150 1200
/usr/bin/time -l $UV wc_copy.py 2,3,5,7 23150 1600 6 2> logs/run8.time | cut -c1-900
rg "real|maximum resident" logs/run8.time
echo RUN8 DONE $(date -u +%H:%M:%S)
