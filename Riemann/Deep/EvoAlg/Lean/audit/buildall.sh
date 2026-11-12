#!/bin/bash
# usage: bash .scratch/evolean/buildall.sh -- clean rebuild of Riemann/Deep/EvoAlg/Lean/*.lean, one module at a time,
# each under the 290 MB watchdog (mk.sh); stops at the first failure; peaks to .scratch/evolean/peaks.log
cd /Users/adriansager/playground/riemann
rm -rf .scratch/evolean/olean
echo "$(date -u +%H:%M) CLEAN REBUILD" >> .scratch/evolean/peaks.log
for m in Defs Basic Lemma2 Diag TwinI Cor41 TheoremA CorA Examples CorA2 CorA3 TheoremN TheoremI2 Lemma2M Check; do
  bash .scratch/evolean/mk.sh $m > /dev/null 2>&1
  rc=$?
  tail -1 .scratch/evolean/peaks.log
  if [ $rc -ne 0 ]; then echo "FAILED at $m"; exit 1; fi
done
echo "ALL OK"
