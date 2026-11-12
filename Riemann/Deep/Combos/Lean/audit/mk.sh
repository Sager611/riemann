#!/bin/bash
# usage: mk.sh Name -- compile Riemann/Deep/Combos/Lean/Name.lean with lean -j1 directly (LEAN_PATH from
# lake env, computed once in .scratch/rclean/leanpath.txt), under .scratch/memcap.sh 290 MB;
# olean to .scratch/rclean/olean; log .scratch/rclean/build_Name.log (peak line from memcap)
cd /Users/adriansager/playground/riemann
O=$PWD/.scratch/rclean/olean
mkdir -p $O/Riemann/Deep/Combos/Lean
n=$1
LOG=.scratch/rclean/build_$n.log
LP=$(cat .scratch/rclean/leanpath.txt)
bash .scratch/memcap.sh 290 1500 env LEAN_PATH="$O:$LP" lean -j1 -o $O/Riemann/Deep/Combos/Lean/$n.olean Riemann/Deep/Combos/Lean/$n.lean > $LOG 2>&1
echo "EXIT=$?" >> $LOG
echo "$(date -u +%H:%M) $n: $(rg 'memcap:' $LOG)" >> .scratch/rclean/peaks.log
head -80 $LOG

