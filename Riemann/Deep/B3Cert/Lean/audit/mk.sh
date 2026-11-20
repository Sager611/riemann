#!/bin/bash
# b3lean build: bash Riemann/Deep/B3Cert/Lean/audit/mk.sh <dir> <Module>  (e.g. Riemann/Deep/B3Cert/Lean Arch)
# lean -j1 directly (no lake build), oleans to .scratch/b3lean/olean, memcap 290 MB / 1800 s, /usr/bin/time -l.
cd /Users/adriansager/playground/riemann
O=/Users/adriansager/playground/riemann/.scratch/b3lean/olean
LP="$O:$(cat Riemann/Deep/B3Cert/Lean/audit/leanpath.txt)"
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
d=$1; m=$2
mkdir -p $O/$d
log=Riemann/Deep/B3Cert/Lean/audit/logs/log_$m.txt
t0=$(date +%s)
LEAN_PATH="$LP" OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 /usr/bin/time -l $LEANBIN -j1 -o $O/$d/$m.olean $d/$m.lean > $log 2>&1
rc=$?
pk=$(rg -o '[0-9]+  peak memory footprint' $log | awk '{printf "%d", $1/1048576}')
mc=$(rg -o 'memcap: peak [0-9]+ MB.*' $log)
echo "$(date -u +%H:%M:%S) $d/$m rc=$rc $(( $(date +%s) - t0 ))s footprint=${pk}MiB $mc" | tee -a Riemann/Deep/B3Cert/Lean/audit/peaks.log
awk '/^ *[0-9.]+ real /{exit} {print}' $log | rg -v '^memcap' | head -${LINES_MAX:-40} | cut -c1-220
exit $rc
