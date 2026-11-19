#!/bin/bash
# explore52: how much of the level-2 gain sits strictly inside the RS range (|lam| <= N-1)? on-line (RH-world) lattice LPs.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for NM in "6 24" "6 36" "7 21" "7 28"; do
  for S in P P3o P3 L4o L4; do
    bash ../../../.scratch/memcap.sh 290 1500 uv run --no-project --with numpy --with scipy python3 npa_cg.py $NM 0.1 $S on 2>&1 | rg 'RESULT|memcap|INFEAS|cap|Error' | sed 's/depths=.*: min/ min/'
  done
done
