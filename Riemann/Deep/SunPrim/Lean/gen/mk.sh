#!/bin/bash
# usage: bash .scratch/sunlean/mk.sh Name  -- compile Riemann/Deep/SunPrim/Lean/Name.lean with lean -j1 (no lake process),
# LEAN_PATH = own olean dir + lake's LEAN_PATH (leanpath.txt), under the 290 MB watchdog; olean to .scratch/sunlean/olean
# probe mode: bash .scratch/sunlean/mk.sh -p path/to/file.lean   (no -o)
cd /Users/adriansager/playground/riemann
O=$PWD/.scratch/sunlean/olean
LP="$O:$(cat .scratch/sunlean/leanpath.txt)"
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
t0=$(date +%s)
if [ "$1" = "-p" ]; then
  f=$2; n=probe_$(basename $f .lean); LOG=.scratch/sunlean/build_$n.log
  LEAN_PATH="$LP" bash .scratch/memcap.sh 290 ${3:-1800} /usr/bin/time -l $LEANBIN -j1 $f > $LOG 2>&1
else
  n=$1; LOG=.scratch/sunlean/build_$n.log
  LEAN_PATH="$LP" bash .scratch/memcap.sh 290 ${2:-1800} /usr/bin/time -l $LEANBIN -j1 -o $O/Riemann/Deep/SunPrim/Lean/$n.olean Riemann/Deep/SunPrim/Lean/$n.lean > $LOG 2>&1
fi
rc=$?
PF=$(rg 'peak memory footprint' $LOG | awk '{printf "%d", $1/1048576}')
echo "$(date -u +%H:%M:%S) $n rc=$rc $(( $(date +%s) - t0 ))s $(rg 'memcap:' $LOG) footprint=${PF}MB" >> .scratch/sunlean/peaks.log
rg -v 'real |user |sys |maximum resident|average shared|average unshared|page reclaims|page faults|swaps|block input|block output|messages sent|messages received|signals received|voluntary context|involuntary context|instructions retired|cycles elapsed|peak memory' $LOG | head -c 6000
echo "footprint ${PF} MB"
echo "mk $n rc=$rc $(( $(date +%s) - t0 ))s"
exit $rc
