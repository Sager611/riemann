#!/bin/bash
# Clean rebuild of all Orders Lean modules, one compile at a time (run from repo root).
cd /Users/adriansager/playground/riemann
export PYTHONDONTWRITEBYTECODE=1
d=Riemann/Deep/Orders/Lean
shasum $d/*.lean > .scratch/orderslean/sha_before.txt
OMP_NUM_THREADS=1 python3 $d/gen.py > /dev/null && OMP_NUM_THREADS=1 python3 $d/gen2.py > /dev/null
shasum $d/*.lean > .scratch/orderslean/sha_after.txt
if cmp -s .scratch/orderslean/sha_before.txt .scratch/orderslean/sha_after.txt; then echo "generators reproduce the data files"; else echo "GENERATED FILES DIFFER"; fi
rm -rf .scratch/orderslean/olean
for m in Defs Ext Bulk Patch Induct D_pxmd D_pxdm D_pdxm D_pdmx D_xpmd D_xpdm D_xdpm D_xdmp D_dxpm D_dpxm D_dpmx Search Fail Small1 Small2 Small3 Small4 Small5 Small6 Main Check; do
  bash .scratch/orderslean/mk.sh $m 1800 > /dev/null
  tail -1 .scratch/orderslean/peaks.log
done
echo DONE

