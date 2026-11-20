#!/bin/bash
# dbnrec2 regen_r4.sh: regenerate the .dat files of tranche 1's row-4 shells 00-05 (prep_end.py, same args) and compare the prep.json line byte for byte.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
ARC=.scratch/dbnrec2/dl/Proof_and_Replay__Lambda_Peer_Review__verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json
AFF=.scratch/dbnrec2/dl/Proof_and_Replay__Lambda_Peer_Review__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
O=Riemann/Deep/DBNRec/logs; L=Riemann/Deep/DBNRec2/logs
us=(0 0.01 0.02 0.035 0.05 0.07 0.09)
for j in 0 1 2 3 4 5; do
  tag=r4_0$j; T=0.834173; [ $j = 0 ] && T=0.832387
  bash .scratch/memcap.sh 290 600 uv run --no-project --with mpmath --with python-flint python3 Riemann/Deep/DBNRec2/code/prep_end.py $ARC $AFF 4 ${us[$j]} ${us[$((j+1))]} $T .scratch/dbnrec2/$tag.dat > .scratch/dbnrec2/$tag.regen.json 2> .scratch/dbnrec2/$tag.regen.err
  if [ "$(head -1 .scratch/dbnrec2/$tag.regen.json)" = "$(head -1 $O/$tag.prep.json)" ]; then echo "$tag prep line IDENTICAL $(tail -1 .scratch/dbnrec2/$tag.regen.err)"; else echo "$tag prep line DIFFERS"; fi
  cp $O/$tag.prep.json $O/$tag.bb.out $O/$tag.unres.txt $L/
done
