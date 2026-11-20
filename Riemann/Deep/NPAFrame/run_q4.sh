#!/bin/bash
# npaframe queue 4: off-line immunity of fixed-theta sets on lattices (Z_M, depths, rays); waits for queue 3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f 'run_q2.sh|run_q3.sh' > /dev/null; do sleep 10; done
D=0.05,0.2,0.6,1.5,4.0
lrun() { tag=L_N$1_M$2_$4_$5${6:+$6}; /usr/bin/time -l bash ../../../.scratch/memcap.sh 290 1500 uv run --no-project --with numpy --with scipy python3 npa_cg_lean.py "$@" > logs/$tag.log 2>&1; echo "$tag :: $(rg 'RESULT' logs/$tag.log | cut -c1-160) :: $(rg -o 'memcap.*' logs/$tag.log)" >> logs/summary.log; }
lrun 6 18 $D P on
lrun 6 18 $D T3 on
lrun 6 18 $D T3 off -rays
lrun 6 18 $D T4 on
lrun 6 18 $D T4 off -rays
lrun 7 21 $D P on
lrun 7 21 $D T5 on
lrun 7 21 $D T5 off -rays
