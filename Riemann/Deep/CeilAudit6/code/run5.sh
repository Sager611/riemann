#!/bin/bash
# ceilaudit6 run5: zero-count convergence of B at the largest N (K doubled).
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy --with mpmath python3"
echo "== zeta K=800"; /usr/bin/time -l $UV b_zeta.py 800 1000,2000,4000 zk800 d2z4000.npy 2>&1 | rg "DELTA|maximum resident|done"
echo "== DH K=600"; /usr/bin/time -l $UV b_dh_k.py 2800 600 d2dh2800.npy 2>&1 | rg "DELTA|maximum resident|done"
echo "== chi K=600"; /usr/bin/time -l $UV b_chi_k.py 3000 900 600 d2c3000.npy 2>&1 | rg "DELTA|maximum resident|done"
echo RUN5 DONE $(date -u +%H:%M:%S)
