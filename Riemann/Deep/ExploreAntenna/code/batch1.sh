#!/bin/bash
# sequential batch of moved-pair thresholds (one process at a time)
cd "$(dirname "$0")"; mkdir -p ../logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
run() { perl -e 'alarm 1500; exec @ARGV' uv run --no-project --with python-flint --with mpmath --with numpy python3 movedpair.py "$@" ; }
run 1 0.03 1.2 6 6 > ../logs/mp_k1_d0.03.log 2>&1
run 1 0.003 1.2 8 6 > ../logs/mp_k1_d0.003.log 2>&1
run 5 0.3 1.5 8 6 > ../logs/mp_k5_d0.3.log 2>&1
run 5 0.03 1.5 9 6 > ../logs/mp_k5_d0.03.log 2>&1
run 10 0.3 2.5 10 6 > ../logs/mp_k10_d0.3.log 2>&1
run 10 0.03 2.5 12 6 > ../logs/mp_k10_d0.03.log 2>&1
run 20 0.3 4 14 6 > ../logs/mp_k20_d0.3.log 2>&1
run 5 0.003 1.5 10 6 > ../logs/mp_k5_d0.003.log 2>&1
run 20 0.03 4 16 6 > ../logs/mp_k20_d0.03.log 2>&1
echo DONE > ../logs/batch1.done
