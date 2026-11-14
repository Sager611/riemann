#!/bin/bash
# certify (K) for all six kernels, one process at a time, T1 = 100 then 2000 if the tail fails (as the paper's driver)
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1
for k in $KERNELS; do
  for T1 in 100 2000; do
    echo "== $k T1=$T1 $(date -u +%H:%M:%S)"
    bash .scratch/memcap.sh 290 900 uv run --no-project --with python-flint python3 Riemann/Deep/STCAudit/certify_K.py $k $T1 800
    rc=$?; echo "rc=$rc"
    [ $rc -eq 0 ] && break
  done
done
echo "== done $(date -u +%H:%M:%S)"
