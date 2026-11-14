#!/bin/bash
# run_all.sh (explore48): rerun of all onset batches after the precision fix (one process at a time, memcap 290 MB).
cd "$(dirname "$0")"
M=../../../.scratch/memcap.sh
U="uv run --no-project --with python-flint python3"
bash $M 290 900 $U t4_onset.py val 3 200 > logs/t4_val.log 2>&1
bash $M 290 1500 $U t4_onset.py h1 3 300 > logs/t4_h1.log 2>&1
bash $M 290 1700 $U t4_onset.py h2 3 150 0 16 > logs/t4_h2.log 2>&1
bash $M 290 1700 $U t4_onset.py rect 3 240 0 7 > logs/t4_rect_0.log 2>&1
bash $M 290 1750 $U t4_onset.py rhomb 3 240 0 13 > logs/t4_rhomb.log 2>&1
echo DONE > logs/run_all.done
