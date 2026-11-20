#!/bin/bash
# ceilaudit6 run2: precision checks (noise on G) + Burnol/Christoffel B(N) for zeta, DH, chi_-4 at large N. One process at a time.
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
KAP=$(cat kap.txt); UV="uv run --no-project --with numpy --with scipy --with mpmath python3"
run() { /usr/bin/time -l bash ../memcap.sh 290 1500 "$@" 2>&1 | rg -v "^ +[0-9]+  (page|block|vol|invol|signals|messages|instructions|cycles|peak|average|swaps)" ; }
echo "== zeta precision"; run ./gram zeta 4000 Gz.bin
run $UV dN_check.py zeta 4000 Gz.bin d2z4000.npy 1e-13
run $UV dN_check.py zeta 4000 Gz.bin d2z4000.npy 1e-11
rm -f Gz.bin
echo "== chi4 precision"; run ./gram per 4 2000 Gc.bin 0 1 1 0
run $UV dN_check.py per 2000 Gc.bin d2c3000.npy 1e-13 4 0 1 1 0
rm -f Gc.bin
echo "== B zeta"; run $UV b_zeta.py 400 10,100,200,400,600,800,1000,1280,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2560,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000 zbig d2z4000.npy
echo "== B DH"; run $UV b_dh.py 2800 300 d2dh2800.npy
echo "== B chi4"; run $UV b_chi.py 3000 900 300 d2c3000.npy
echo RUN2 DONE $(date -u +%H:%M:%S)
