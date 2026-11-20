#!/bin/bash
# ceilaudit6 run1: big NB Grams (zeta 4000, chi_-4 3000, DH 2800), one process at a time, memcap 290 MB.
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
KAP=$(cat kap.txt)
run() { /usr/bin/time -l bash ../memcap.sh 290 1200 "$@" 2>&1 | rg -v "^ +[0-9]+  (page|block|vol|invol|signals|messages|instructions|cycles|peak|average|swaps)" ; }
run ./gram zeta 4000 Gz4000.bin
run uv run --no-project --with numpy --with scipy python3 dN.py zeta 4000 Gz4000.bin d2z4000.npy
rm -f Gz4000.bin
run ./gram per 4 3000 Gc3000.bin 0 1 1 0
run uv run --no-project --with numpy --with scipy python3 dN.py per 3000 Gc3000.bin d2c3000.npy 4 0 1 1 0
rm -f Gc3000.bin
run ./gram per 5 2800 Gdh2800.bin 0 1 $KAP 1 0
run uv run --no-project --with numpy --with scipy python3 dN.py per 2800 Gdh2800.bin d2dh2800.npy 5 0 1 $KAP 1 0
rm -f Gdh2800.bin
echo RUN1 DONE $(date -u +%H:%M:%S)
