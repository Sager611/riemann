#!/bin/bash
# dpcblean build (copy of dpclean's mk.sh): bash .scratch/dpcblean/mk.sh <Module> [dir]   (default dir Riemann/Deep/DPCBeyond/Lean)
# lean -j1 directly (no lake), oleans to .scratch/dpcblean/olean, under memcap 290 MB / 1800 s, /usr/bin/time -l.
cd /Users/adriansager/playground/riemann
O=/Users/adriansager/playground/riemann/.scratch/dpcblean/olean
LP="$O:$(cat .scratch/dpcblean/leanpath.txt)"
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
m=$1; d=${2:-Riemann/Deep/DPCBeyond/Lean}
mkdir -p $O/$d
log=.scratch/dpcblean/log_$m.txt
t0=$(date +%s)
LEAN_PATH="$LP" OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 /usr/bin/time -l $LEANBIN -j1 -o $O/$d/$m.olean $d/$m.lean > $log 2>&1
rc=$?
pk=$(rg -o '[0-9]+  peak memory footprint' $log | awk '{printf "%d", $1/1048576}')
mc=$(rg -o 'memcap: peak [0-9]+ MB.*' $log)
echo "$(date -u +%H:%M:%S) $m rc=$rc $(( $(date +%s) - t0 ))s footprint=${pk}MiB $mc" | tee -a .scratch/dpcblean/peaks.log
awk '/^ *[0-9.]+ real /{exit} {print}' $log | rg -v '^memcap' | head -${LINES_MAX:-80}
exit $rc
