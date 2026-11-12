#!/bin/bash
# sunlean clean rebuild of all compiled modules, in dependency order, one at a time (memcap 290 MB).
cd /Users/adriansager/playground/riemann
rm -rf .scratch/sunlean/olean
for m in Graph Arith Prop4 Defs Sieve CaseA1 TheoremA Small Cycles Cycles9 CheckA CheckS FFCoreC; do
  LINES_MAX=40 bash .scratch/sunlean/mk.sh $m || echo "FAILED $m"
done
echo DONE
