#!/bin/bash
cd /Users/adriansager/playground/riemann/Riemann/Deep/ExploreSchurWeyl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1
MC=../../../.scratch/memcap.sh
for a in "12 96 tri SCS" "12 96 pair SCS" "16 128 tri SCS" "16 128 pair SCS"; do
  set -- $a
  bash $MC 290 1200 uv run --no-project --with numpy --with scipy --with cvxpy --with clarabel python3 tri_sdp.py $1 $2 $3 $4 2>&1 | rg -v 'Warning|warn' | tail -2 >> logs/tri_sdp.log
done
echo QUEUE2_DONE >> logs/tri_sdp.log
