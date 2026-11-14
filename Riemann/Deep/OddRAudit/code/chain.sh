#!/bin/bash
# sequential reruns, one process at a time, each under memcap 290 MB
cd /Users/adriansager/playground/riemann/.scratch/oddraudit/run
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY="uv run --no-project --with python-flint --with mpmath --with numpy python3"
for c in "gap.py" "dfs_gap.py" "oddca.py 1e6" "oddca_seg.py 1e6" "oddca.py 5e7" "oddca_seg.py 5e7" "transfer.py 5e8 1.43e47 0.87338196661320752" "transfer.py 5e8 1.42e47 7/8" "transfer.py 5e8 6.46e46 1" "oddca_seg.py 2e8" "oddca_seg.py 5e8"; do
  echo "== $c  ($(date -u +%H:%M:%S))"
  bash ../../memcap.sh 290 1800 /usr/bin/time -l $PY $c > out.tmp 2> err.tmp; rc=$?
  cat out.tmp | tail -12
  rg "memcap|maximum resident|peak memory footprint|real" err.tmp
  echo "rc=$rc"
done
echo "CHAIN DONE $(date -u +%H:%M:%S)"
