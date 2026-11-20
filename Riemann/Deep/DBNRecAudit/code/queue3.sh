#!/bin/bash
# dbnrecaudit queue3.sh: after queue2.sh, regenerate 8 shells' .dat (shell 0 of rows 1-7, DBNRec2 r4_09 = largest K) with the repo's prep_end.py
# into /tmp/dbnra_keep (prep.json line must equal the logged one), then phase_arb.py on these + r8_00.dat. One process at a time, memcap 290 MB.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code; K=/tmp/dbnra_keep
UV="uv run --no-project --with mpmath --with numpy --with python-flint python3"
DL=.scratch/dbnrecaudit/dl/Proof_and_Replay__Lambda_Peer_Review__verification
ARC=${DL}__arcs__ACTUAL_67_ROW_ARC_INPUTS.json; AFF=${DL}__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
while pgrep -f 'DBNRecAudit/code/queue2?.sh' > /dev/null; do sleep 5; done
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/queue.log; }
specs="$K/r8_00.dat:Riemann/Deep/DBNRec2/logs/r8_00.prep.json"
for dt in DBNRec:r1_00 DBNRec:r2_00 DBNRec:r3_00 DBNRec:r4_00 DBNRec2:r5_00 DBNRec2:r6_00 DBNRec2:r7_00 DBNRec2:r4_09; do
  d=${dt%%:*}; t=${dt##*:}; pj=Riemann/Deep/$d/logs/$t.prep.json
  args=$(python3 -c "import json;from fractions import Fraction as F;p=json.loads(open('$pj').readline());print(p['row'],'%.9f'%float(F(p['u0'])),'%.9f'%float(F(p['u1'])),'%.9f'%float(F(p['T'])))")
  log "regen $d/$t args $args"
  bash .scratch/memcap.sh 290 600 $UV Riemann/Deep/DBNRec2/code/prep_end.py $ARC $AFF $args $K/${d}_$t.dat > $K/${d}_$t.json 2> $K/${d}_$t.err
  if head -1 $K/${d}_$t.json | cmp -s - <(head -1 $pj); then log "regen $d/$t prep.json identical"; specs="$specs $K/${d}_$t.dat:$pj"; else log "regen $d/$t prep.json DIFFERS"; fi
done
log "start phase_arb"
/usr/bin/time -l bash .scratch/memcap.sh 290 900 $UV $C/phase_arb.py $specs > $L/phase_arb.out 2> $L/phase_arb.err
log "end phase_arb rc=$? $(rg -o 'memcap: peak [0-9]+ MB.*' $L/phase_arb.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/phase_arb.err | tail -1)"
rm -f $K/DBNRec*_r*.dat
log "queue3 done"
