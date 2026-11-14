#!/bin/bash
# e1lean clean rebuild of all modules, one at a time (memcap 290 MB / 1800 s each).
cd /Users/adriansager/playground/riemann
rm -rf .scratch/e1lean/olean && mkdir -p .scratch/e1lean/olean
: > .scratch/e1lean/verify.log
for m in Defs Detector AlgebraH AlgebraC Numeric Kappa Arch Parity Bounds2 Check; do
  LINES_MAX=1 bash .scratch/e1lean/mk.sh $m >> .scratch/e1lean/verify.log 2>&1 || { echo "FAIL $m" >> .scratch/e1lean/verify.log; exit 1; }
done
echo "ALL OK" >> .scratch/e1lean/verify.log
