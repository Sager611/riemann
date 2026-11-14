#!/bin/bash
# Reruns every script of Riemann/Deep/B2Below in sequence (one process at a time) under the memory watchdog; output -> run.log
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
UV="uv run --no-project --with mpmath --with numpy --with sympy python3"
{
for s in chain.py system.py infimum.py chain_limit.py ceiling.py levers.py perturb.py ode.py check_identity.py poly143.py theta_inf.py obstruction.py; do
  echo "=== $s ($(date -u +%H:%M:%S) UTC)"
  if [ $s = chain.py ] || [ $s = chain_limit.py ] || [ $s = check_identity.py ]; then cmd="python3 $s"; else cmd="$UV $s"; fi
  bash ../../../.scratch/memcap.sh 290 900 /usr/bin/time -l $cmd 2>&1 | rg -v '^\s+[0-9]+ +(involuntary|voluntary|signals|messages|block|page|swaps|average|instructions|cycles)'
done
} > run.log 2>&1
