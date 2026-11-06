#!/bin/bash
cd "$(dirname "$0")"; mkdir -p ../logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
U="uv run --no-project --with python-flint --with mpmath --with numpy python3"
perl -e 'alarm 1200; exec @ARGV' $U dh.py 42 44 46 > ../logs/dh_second.log 2>&1
perl -e 'alarm 1500; exec @ARGV' $U movedpair.py 30 0.3 5 12 6 > ../logs/mp_k30_d0.3.log 2>&1
NMUL=1.6 perl -e 'alarm 1500; exec @ARGV' $U movedpair.py 10 0.3 3.5 5.5 5 > ../logs/mp_k10_d0.3_N16.log 2>&1
perl -e 'alarm 1800; exec @ARGV' $U movedpair.py 50 0.3 7 15 6 > ../logs/mp_k50_d0.3.log 2>&1
echo DONE > ../logs/batch2.done
