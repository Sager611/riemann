#!/bin/bash
cd /Users/adriansager/playground/riemann/Riemann/Deep/ExploreSchurWeyl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1
while pgrep -f queue2.sh > /dev/null; do sleep 5; done
MC=../../../.scratch/memcap.sh
for a in "4 48 trig SCS" "4 48 pairg SCS" "8 64 trig SCS" "8 64 pairg SCS"; do
  set -- $a
  bash $MC 290 900 uv run --no-project --with numpy --with scipy --with cvxpy --with clarabel python3 tri_sdp.py $1 $2 $3 $4 2>&1 | rg -v 'Warning|warn' | tail -2 >> logs/tri_sdp.log
done
echo QUEUE3_DONE >> logs/tri_sdp.log
