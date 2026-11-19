#!/bin/bash
# explore52: N=6 column-generation runs (one process at a time, memcap 290 MB, single-threaded); -rays = infinitely deep pairs
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
D=0.03,0.1,0.3,0.6,1.0,2.0,3.5,5.0
for SW in "P off" "P3 on" "P3 off" "L4 on" "L4 off"; do
  set -- $SW
  bash ../../../.scratch/memcap.sh 290 1500 uv run --no-project --with numpy --with scipy python3 npa_cg.py 6 24 $D $1 $2 -rays 2>&1 | rg 'it=|RESULT|memcap|INFEAS|cap|rays|Error' | sed 's/depths=.*: min/ min/'
done
