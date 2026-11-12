#!/bin/bash
cd /Users/adriansager/playground/riemann/Riemann/Deep/GapCeil
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for k in 3.305 3.35 3.50 4.00 5.00; do
  perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l uv run --no-project --with python-flint --with numpy --with scipy python3 runcert_r.py $k 72 2 0.8 1.7 > logs/certR_r2_$k.txt 2>&1
  rg -n 'RESULT|maximum resident' logs/certR_r2_$k.txt | head -2 >> logs/summaryR.txt
done
echo BATCHR DONE >> logs/summaryR.txt
