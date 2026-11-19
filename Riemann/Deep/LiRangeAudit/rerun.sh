#!/bin/bash
# rerun.sh -- LiRangeAudit own reruns, one process at a time, each under memcap 290 MB
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=Riemann/Deep/LiRangeAudit/logs
bash .scratch/memcap.sh 290 900 /usr/bin/time -l uv run --no-project --with mpmath python3 Riemann/Deep/LiRangeAudit/indep_g.py > $L/indep_g.log 2>&1
bash .scratch/memcap.sh 290 1200 /usr/bin/time -l uv run --no-project --with python-flint python3 Riemann/Deep/HP141/Li200k/verify.py > $L/rerun_verify141.log 2>&1
bash .scratch/memcap.sh 290 600 /usr/bin/time -l uv run --no-project --with python-flint --with numpy python3 Riemann/Deep/DirLi/certify.py small > $L/rerun_dirli_small.log 2>&1
bash .scratch/memcap.sh 290 900 /usr/bin/time -l uv run --no-project --with python-flint --with numpy python3 Riemann/Deep/DirLi/certify.py blocks 1.02 301948 307985 > $L/rerun_dirli_block.log 2>&1
echo DONE >> $L/rerun_dirli_block.log
