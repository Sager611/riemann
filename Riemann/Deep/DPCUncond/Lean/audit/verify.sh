#!/bin/bash
# dpclean final verification: every module in dependency order, one at a time, lean -j1 directly,
# memcap 290 MB / 1800 s, /usr/bin/time -l; pplean's FinSum, Model, Liminf first (imported, unmodified).
cd /Users/adriansager/playground/riemann
out=.scratch/dpclean/verify.log; : > $out
for m in FinSum Model Liminf; do
  LINES_MAX=5 bash .scratch/dpclean/mk.sh $m Riemann/Deep/ExplorePointProc/Lean >> $out 2>&1
done
for m in Params Gaps Tail Phi Step7 TheoremU DLam Hairpin Fine Check; do
  LINES_MAX=200 bash .scratch/dpclean/mk.sh $m >> $out 2>&1
done
cat $out

