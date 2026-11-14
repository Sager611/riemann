#!/bin/bash
# run2.sh (explore48): refined reruns (local N0 refinement), off-line birth distribution, large y. One process at a time.
cd "$(dirname "$0")"
M=../../../.scratch/memcap.sh
U="uv run --no-project --with python-flint python3"
while [ ! -f logs/run_all.done ]; do sleep 5; done
bash $M 290 1200 $U t4_onset.py h1 3 400 > logs/t4_h1r.log 2>&1
bash $M 290 1200 $U t4_onset.py redo 3 300 > logs/t4_redo.log 2>&1
bash $M 290 900 $U t6_dist.py rect10 5 > logs/t6_rect10.log 2>&1
bash $M 290 900 $U t6_dist.py D427 5 > logs/t6_D427.log 2>&1
bash $M 290 900 $U t6_dist.py D163 4 > logs/t6_D163.log 2>&1
bash $M 290 1800 $U t4_onset.py big 4 1500 > logs/t4_big.log 2>&1
echo DONE > logs/run2.done
