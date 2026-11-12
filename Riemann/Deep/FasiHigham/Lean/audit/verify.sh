#!/bin/bash
# Clean rebuild of all FasiHigham Lean modules, one at a time, under memcap 290 MB (run from repo root).
set -e
rm -rf .scratch/fasilean/olean
d=Riemann/Deep/FasiHigham/Lean
for m in Defs Core KeyIneq Bound Root PS Main MainPS Quant Quant2; do
  bash .scratch/fasilean/mk.sh $d/$m.lean .scratch/fasilean/olean/$d/$m.olean > /dev/null
  tail -1 .scratch/fasilean/peaks.log
  if rg -q 'error' .scratch/fasilean/log_$m.txt; then echo "ERROR in $m"; exit 1; fi
done
bash .scratch/fasilean/mk.sh $d/Check.lean > .scratch/fasilean/check_out.txt
tail -1 .scratch/fasilean/peaks.log
rg -c 'depends on axioms: \[propext, Classical.choice, Quot.sound\]' .scratch/fasilean/check_out.txt
