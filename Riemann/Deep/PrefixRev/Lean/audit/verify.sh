#!/bin/bash
# Clean rebuild of Riemann/Deep/PrefixRev/Lean in dependency order, one module at a time, under memcap 290 MB.
cd /Users/adriansager/playground/riemann
rm -rf .scratch/prefixlean/olean
for m in Defs CycleA CycleB CycleC Parity IsCycleP Small BlocksA BlocksB Main PrimA PrimB1 PrimB2 Jones Orbits CheckA CheckB; do
  bash .scratch/prefixlean/mk.sh $m 1800 3
done
echo "axiom lines (std): $(rg -c 'depends on axioms: \[propext, Classical.choice, Quot.sound\]|depends on axioms: \[propext, Quot.sound\]' .scratch/prefixlean/logs/CheckA.log .scratch/prefixlean/logs/CheckB.log | tr '\n' ' ')"
echo "nonstd axiom lines: $(rg 'depends on axioms' .scratch/prefixlean/logs/CheckA.log .scratch/prefixlean/logs/CheckB.log | rg -v 'propext, Classical.choice, Quot.sound\]|propext, Quot.sound\]' | wc -l)"
echo "errors: $(rg -c 'error' .scratch/prefixlean/logs/*.log | rg -v ':0$' | tr '\n' ' ')"
echo VERIFY_DONE
