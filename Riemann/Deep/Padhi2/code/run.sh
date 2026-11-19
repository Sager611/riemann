#!/bin/bash
# usage: bash code/run.sh <timeout_s> <name> <script.py> [args]  -> logs/<name>.log, logs/<name>.err (memcap peak + /usr/bin/time -l)
cd "$(dirname "$0")/.."
to=$1; name=$2; shift 2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/usr/bin/time -l bash ../../../.scratch/memcap.sh 290 "$to" uv run --no-project --with python-flint --with mpmath --with numpy --with scipy python3 "$@" > logs/$name.log 2> logs/$name.err
rc=$?
echo "rc=$rc"; rg "memcap|maximum resident|peak memory footprint" logs/$name.err | head -3
exit $rc
