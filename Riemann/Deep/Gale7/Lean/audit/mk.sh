#!/bin/bash
# usage: bash Riemann/Deep/Gale7/Lean/audit/mk.sh Module -- compile Riemann/Deep/Gale7/Lean/Module.lean with lean -j1
# (LEAN_PATH = own olean dir + lake's LEAN_PATH), under the 290 MB / 1800 s watchdog; oleans to .scratch/galelean/olean
cd /Users/adriansager/playground/riemann
O=$PWD/.scratch/galelean/olean
mkdir -p $O/Riemann/Deep/Gale7/Lean .scratch/galelean
[ -s .scratch/galelean/leanpath.txt ] || lake env printenv LEAN_PATH > .scratch/galelean/leanpath.txt
m=$1
LOG=.scratch/galelean/build_$m.log
LP=$(cat .scratch/galelean/leanpath.txt)
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
t0=$(date +%s)
LEAN_PATH="$O:$LP" bash .scratch/memcap.sh 290 1800 /usr/bin/time -l $LEANBIN -j1 -o $O/Riemann/Deep/Gale7/Lean/$m.olean Riemann/Deep/Gale7/Lean/$m.lean > $LOG 2>&1
rc=$?
fp=$(rg -o "[0-9]+  peak memory footprint" $LOG | awk '{printf "%.0f", $1/1048576}')
echo "$(date -u +%H:%M) $m rc=$rc $(( $(date +%s) - t0 ))s footprint=${fp}MB $(rg -o 'memcap: peak [0-9]+ MB' $LOG)"
rg -v "deprecated|push_neg|^$" $LOG | rg "error|depends on axioms" | head -40
exit $rc
