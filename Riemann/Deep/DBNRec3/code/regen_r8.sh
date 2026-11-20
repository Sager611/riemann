#!/bin/bash
# dbnrec3 regen_r8.sh: regenerate the .dat files of tranche 2's row-8 shells 00-05 (prep_end.py, same args as tranche 2's dorow), compare the prep.json line
# byte for byte with DBNRec2/logs (read-only), and copy those shells' prep.json/bb.out/unres.txt into DBNRec3/logs (finished inputs).
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
ARC=.scratch/dbnrec3/dl/Proof_and_Replay__Lambda_Peer_Review__verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json
AFF=.scratch/dbnrec3/dl/Proof_and_Replay__Lambda_Peer_Review__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
O=Riemann/Deep/DBNRec2/logs; L=Riemann/Deep/DBNRec3/logs
us=(0 0.01 0.02 0.035 0.05 0.07 0.09)
for j in 0 1 2 3 4 5; do
  tag=r8_0$j; T=0.913695; [ $j = 0 ] && T=0.913118
  bash .scratch/memcap.sh 290 600 uv run --no-project --with mpmath --with python-flint python3 Riemann/Deep/DBNRec3/code/prep_end.py $ARC $AFF 8 ${us[$j]} ${us[$((j+1))]} $T .scratch/dbnrec3/$tag.dat > .scratch/dbnrec3/$tag.regen.json 2> .scratch/dbnrec3/$tag.regen.err
  if [ "$(head -1 .scratch/dbnrec3/$tag.regen.json)" = "$(head -1 $O/$tag.prep.json)" ]; then echo "$tag prep line IDENTICAL $(tail -1 .scratch/dbnrec3/$tag.regen.err)"; else echo "$tag prep line DIFFERS"; fi
  cp $O/$tag.prep.json $O/$tag.bb.out $O/$tag.unres.txt $L/
done
