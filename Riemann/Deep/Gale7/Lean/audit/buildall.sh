#!/bin/bash
# usage: bash Riemann/Deep/Gale7/Lean/audit/buildall.sh -- clean rebuild, one module at a time, each under the watchdog
cd /Users/adriansager/playground/riemann
rm -rf .scratch/galelean/olean
for m in Lemma Desc Fam Defs Crit Simplex3 Simplex4 Q7 Q8 N7 N8 Main Check; do
  bash Riemann/Deep/Gale7/Lean/audit/mk.sh $m || { echo "FAILED at $m"; exit 1; }
done
echo "ALL OK"
