#!/bin/bash
# explore52 queue (waits for run_inside.sh): (1) N = 6 open-top sets, off-line world; (2) closed level-2 block L4c (NPA principal
# block: every off-diagonal word with both diagonal partners) on/off + rays at N = 4, 5; (3) continuous-position colgen N = 8.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while pgrep -f run_inside.sh > /dev/null; do sleep 5; done
for S in P3o L4o; do
  bash ../../../.scratch/memcap.sh 290 1200 uv run --no-project --with numpy --with scipy python3 npa_cg.py 6 24 0.03,0.1,0.3,0.6,1.0,2.0,3.5,5.0 $S off 2>&1 | rg 'RESULT|memcap|INFEAS|cap|Error' | sed 's/depths=.*: min/ min/' >> logs/cg_N6_M24_opentop_off.log
done
for NM in "4 48" "5 40"; do for W in on off; do
  bash ../../../.scratch/memcap.sh 290 1200 uv run --no-project --with numpy --with scipy python3 npa_cg.py $NM 0.01,0.03,0.1,0.2,0.35,0.6,1.0,1.5,2.2,3.0,4.0,5.0 L4c $W -rays 2>&1 | rg 'RESULT|memcap|INFEAS|cap|Error' | sed 's/depths=.*: min/ min/' >> logs/cg_L4c.log
done; done
for S in P P3o; do
  bash ../../../.scratch/memcap.sh 290 1100 uv run --no-project --with numpy --with scipy python3 colgen_open.py 8 $S 60 384 1 > logs/colgen_N8_$S.log 2>&1
done
