#!/bin/bash
# usage: .scratch/oddrlean/mk.sh Module  -- compile Riemann/Deep/OddRobin/Lean/Module.lean with lean directly (no lake process),
# LEAN_PATH = private olean dir + lake's LEAN_PATH (saved once in leanpath.txt); under the 290 MB watchdog; peaks to peaks.log
cd /Users/adriansager/playground/riemann
m=$1
O=/Users/adriansager/playground/riemann/.scratch/oddrlean/olean
mkdir -p $O/Riemann/Deep/OddRobin/Lean
LP="$O:$(cat .scratch/oddrlean/leanpath.txt)"
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
t0=$(date +%s)
LEAN_PATH="$LP" bash .scratch/memcap.sh 290 1800 $LEANBIN -j1 -o $O/Riemann/Deep/OddRobin/Lean/$m.olean Riemann/Deep/OddRobin/Lean/$m.lean > .scratch/oddrlean/build_$m.log 2>&1
rc=$?
echo "$(date -u +%H:%M) $m rc=$rc $(( $(date +%s) - t0 ))s $(rg 'memcap:' .scratch/oddrlean/build_$m.log)" >> .scratch/oddrlean/peaks.log
cat .scratch/oddrlean/build_$m.log | head -60
echo "mk $m rc=$rc $(( $(date +%s) - t0 ))s"
exit $rc

