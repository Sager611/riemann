#!/bin/bash
# Sequential runs (one process at a time), each under the 290 MB / 900 s watchdog.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
MC=../../../../.scratch/memcap.sh
UV="uv run --no-project --with mpmath --with numpy --with scipy --with sympy python3"
for m in "$@"; do
  case $m in
    orig30) BFIX_MODE=x bash $MC 290 900 $UV harness.py constants_orig.py tables.json ../logs/orig30.json 30 > ../logs/orig30.log 2> ../logs/orig30.err ;;
    *30) mm=${m%30}; BFIX_MODE=$mm bash $MC 290 900 $UV harness.py constants_fix.py tables.json ../logs/$m.json 30 > ../logs/$m.log 2> ../logs/$m.err ;;
    c5only16) BFIX_G=23.5604 BFIX_MODE=c5only bash $MC 290 900 $UV harness.py constants_fix.py tables.json ../logs/c5only16.json > ../logs/c5only16.log 2> ../logs/c5only16.err ;;
    fix16) BFIX_G=23.5604 BFIX_MODE=fix bash $MC 290 900 $UV harness.py constants_fix.py tables.json ../logs/fix16.json > ../logs/fix16.log 2> ../logs/fix16.err ;;
    *) BFIX_MODE=$m bash $MC 290 900 $UV harness.py constants_fix.py tables.json ../logs/$m.json > ../logs/$m.log 2> ../logs/$m.err ;;
  esac
  echo "$m rc=$? $(rg -o 'memcap.*' ../logs/$m.err)"
done
