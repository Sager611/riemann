#!/bin/bash
R=/Users/adriansager/playground/riemann
MC="bash $R/.scratch/memcap.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
run() { local name=$1 lim=$2; shift 2; echo "== $name start $(date -u +%H:%M:%S)"; /usr/bin/time -l $MC 290 $lim "$@" > /tmp/stc2/logs/$name.log 2>&1; echo "rc=$? end $(date -u +%H:%M:%S)"; rg 'memcap:' /tmp/stc2/logs/$name.log; }
cd /tmp/stc2/Ingham3 && run Ingham3 1400 uv run --no-project --with python-flint --with mpmath --with sympy python3 cert.py
cd /tmp/stc2/RobinCH && run RobinCH_record 1200 uv run --no-project --with python-flint --with mpmath --with numpy python3 /tmp/stc2/robinch_record.py
cd /tmp/stc2/RobinHeight && run RobinHeight_record 1200 uv run --no-project --with python-flint --with mpmath --with numpy python3 /tmp/stc2/robinheight_record.py
echo "== all done $(date -u +%H:%M:%S)"
