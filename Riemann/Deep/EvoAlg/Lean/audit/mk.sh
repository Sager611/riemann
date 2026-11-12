#!/bin/bash
# usage: bash .scratch/evolean/mk.sh Module -- compile Riemann/Deep/EvoAlg/Lean/Module.lean with lean -j1 directly
# (LEAN_PATH = own olean dir + lake's LEAN_PATH, saved in leanpath.txt), under memcap 290 MB; olean to .scratch/evolean/olean
cd /Users/adriansager/playground/riemann
O=$PWD/.scratch/evolean/olean
mkdir -p $O/Riemann/Deep/EvoAlg/Lean
m=$1
LOG=.scratch/evolean/build_$m.log
LP=$(cat .scratch/evolean/leanpath.txt)
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
t0=$(date +%s)
LEAN_PATH="$O:$LP" bash .scratch/memcap.sh 290 1800 $LEANBIN -j1 -o $O/Riemann/Deep/EvoAlg/Lean/$m.olean Riemann/Deep/EvoAlg/Lean/$m.lean > $LOG 2>&1
rc=$?
echo "rc=$rc" >> $LOG
echo "$(date -u +%H:%M) $m rc=$rc $(( $(date +%s) - t0 ))s $(rg 'memcap:' $LOG)" >> .scratch/evolean/peaks.log
head -100 $LOG
echo "mk $m rc=$rc $(( $(date +%s) - t0 ))s"
exit $rc
