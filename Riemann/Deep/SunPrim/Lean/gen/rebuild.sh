#!/bin/bash
# clean rebuild of all modules in dependency order; summary to .scratch/sunlean/rebuild.log
cd /Users/adriansager/playground/riemann
rm -rf .scratch/sunlean/olean; mkdir -p .scratch/sunlean/olean/Riemann/Deep/SunPrim/Lean
L=.scratch/sunlean/rebuild.log; : > $L
MODS="Groups ListAux Conj310 A228728 FieldsBase Conj49 Conj49b Conj49c MO316143 WangSun Roots RealInst CircAux F8Defs F8Add F8Mul F8DistL F8DistR F8 F8Search F8Leaf0 F8Leaf1 F8Leaf2 F8Leaf3 F8Top F9Defs F9Add F9Mul F9DistL F9DistR F9 F9Search"
for i in $(seq 0 17); do MODS="$MODS F9Leaf$i"; done
MODS="$MODS F9Top A229827"
for m in $MODS; do
  out=$(bash .scratch/sunlean/mk.sh $m 900 2>&1)
  rc=$(echo "$out" | tail -1 | awk '{print $3}')
  fp=$(echo "$out" | rg '^footprint' | awk '{print $2}')
  err=$(rg -c 'error' .scratch/sunlean/build_$m.log)
  warn=$(rg -c 'warning' .scratch/sunlean/build_$m.log)
  ax=$(rg "depends on axioms" .scratch/sunlean/build_$m.log | tr '\n' ' ')
  echo "$m $rc footprint=${fp}MB errors=${err:-0} warnings=${warn:-0} $ax" >> $L
done
echo DONE >> $L
