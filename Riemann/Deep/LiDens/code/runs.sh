#!/bin/bash
# runs.sh -- the certificate runs of REPORT.md, one process at a time (from the repo root)
cd "$(dirname "$0")/../../../.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R="uv run --no-project --with python-flint --with mpmath --with sympy python3 Riemann/Deep/LiDens/code/cert.py"
W=""; [ -f .scratch/memcap.sh ] && W="bash .scratch/memcap.sh 290 1500"     # repo watchdog (290 MB, 25 min) when present
run() { tag=$1; shift; /usr/bin/time -l $W $R "$@" > Riemann/Deep/LiDens/logs/cert_$tag.log 2>&1; }
[ "$1" = "all" ] && run B141 --B 141 --sstar 0.83 --J 176 --falsify
[ "$1" = "best" ] || run B144 --B 144 --sstar 0.83 --J 176
[ "$1" = "best" ] || run B141_plain --B 141 --sstar 0.845 --J 189 --plain --nomin
[ "$1" = "best" ] || run B141_delta1 --B 141 --sstar 0.83 --J 176 --delta1
# best run (non-uniform slabs, finer t-grid; also part 1 and the real-scale falsification)
[ "$1" = "best" ] || [ "$1" = "all" ] && run B141_best --B 141 --mesh 0.80:0.002,0.84:0.0001 --fine --falsify --small
[ "$1" = "best" ] || [ "$1" = "all" ] && run B144_best --B 144 --mesh 0.80:0.002,0.84:0.0001 --fine
