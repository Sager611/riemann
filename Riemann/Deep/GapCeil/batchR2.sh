#!/bin/bash
cd /Users/adriansager/playground/riemann/Riemann/Deep/GapCeil
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while pgrep -f batchR.sh > /dev/null; do sleep 3; done
perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l uv run --no-project --with python-flint --with numpy --with scipy python3 runcert_r.py 3.305 72 2 0.8 1.7 0.4 > logs/certR_r2_3.305.txt 2>&1
rg -n 'RESULT|maximum resident' logs/certR_r2_3.305.txt | head -2 >> logs/summaryR.txt
perl -e 'alarm 600; exec @ARGV' /usr/bin/time -l uv run --no-project --with python-flint python3 verify_q.py > logs/verify_q.txt 2>&1
echo BATCHR2 DONE >> logs/summaryR.txt
