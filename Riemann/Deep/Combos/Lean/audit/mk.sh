#!/bin/bash
# usage: bash .scratch/rclean2/mk.sh Name [file]  -- compile Riemann/Deep/Combos/Lean/Name.lean (or [file], then no olean)
# lean -j1 direct binary, LEAN_PATH = own olean dir + lake env path (leanpath.txt); under memcap 290 MB / 1800 s;
# /usr/bin/time -l gives the exact peak footprint; one line per run in .scratch/rclean2/peaks.log
cd /Users/adriansager/playground/riemann
S=.scratch/rclean2
O=$PWD/$S/olean
n=$1; f=${2:-Riemann/Deep/Combos/Lean/$n.lean}
LOG=$S/build_$n.log
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
if [ -z "$2" ]; then OUT="-o $O/Riemann/Deep/Combos/Lean/$n.olean"; else OUT=""; fi
t0=$(date +%s)
LEAN_PATH="$O:$(cat $S/leanpath.txt)" bash .scratch/memcap.sh 290 1800 /usr/bin/time -l $LEANBIN -j1 $LEANOPTS $OUT $f > $LOG 2>&1
rc=$?
pf=$(rg 'peak memory footprint' $LOG | awk '{printf "%d", $1/1048576}')
echo "$(date -u +%H:%M:%S) $n rc=$rc $(( $(date +%s) - t0 ))s exact=${pf}MB $(rg -o 'memcap:.*' $LOG)" >> $S/peaks.log
rg -v '^\s+[0-9]+\s+[a-z]|real.*user.*sys' $LOG | head -40
exit $rc
