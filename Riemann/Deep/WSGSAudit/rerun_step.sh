#!/bin/sh
# rerun_step.sh c eps NG bits K2: one chain step's two certificates at c (defl: U(c) and lambda_2^even >= eps; odd: lambda_1^odd >= eps)
# from the copy /tmp/wsgs_copy (never the repo folder), fresh scan cache, single-threaded, each under memcap 290 MB.
c=$1; eps=$2; NG=$3; bits=$4; K2=$5
mkdir -p /tmp/wsgs_scan
for kind in defl odd; do
  OMP_NUM_THREADS=1 WS_NG=$NG WS_M=1 WZ_RSCAN=/tmp/ws_rscan WZ_SCAN=/tmp/wsgs_scan bash /Users/adriansager/playground/riemann/.scratch/memcap.sh 290 1200 \
    /usr/bin/time -l uv run --no-project --with mpmath --with numpy --with python-flint python3 /tmp/wsgs_copy/cert_simple.py $c $kind $eps $bits $K2 2>/tmp/wsgs_err_$kind.txt
  rg -o "memcap.*|[0-9]+ +maximum resident set size" /tmp/wsgs_err_$kind.txt | tr '\n' ' '; echo
done
