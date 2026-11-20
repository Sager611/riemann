#!/bin/bash
# dbnrecaudit queue5.sh (17:08 UTC; replaces the rest of queue.sh-queue4.sh, cut for time): wait for the running dense row 2, then dense rows 1, 8
# (cells 10^6 points; wide windows 2*10^5), the Arb phase/domain check (queue3.sh's body), dense row 4, the 32-offset scan of row 8.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code; K=/tmp/dbnra_keep
UV="uv run --no-project --with mpmath --with numpy --with python-flint python3"
DL=.scratch/dbnrecaudit/dl/Proof_and_Replay__Lambda_Peer_Review__verification
ARC=${DL}__arcs__ACTUAL_67_ROW_ARC_INPUTS.json; AFF=${DL}__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/queue.log; }
while pgrep -f 'falsify_dense.py 2' > /dev/null; do sleep 3; done
log "end dense_r2 (queue.sh stopped at 17:08 for time) $(rg -o 'memcap: peak [0-9]+ MB.*' $L/dense_r2.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/dense_r2.err | tail -1)"
step() { name=$1; to=$2; shift 2; log "start $name"
  /usr/bin/time -l bash .scratch/memcap.sh 290 $to "$@" > $L/$name.out 2> $L/$name.err; rc=$?
  log "end $name rc=$rc $(rg -o 'memcap: peak [0-9]+ MB.*' $L/$name.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/$name.err | tail -1)"; }
step dense_r1 1500 $UV $C/falsify_dense.py 1
step dense_r8 1500 $UV $C/falsify_dense.py 8
specs="$K/r8_00.dat:Riemann/Deep/DBNRec2/logs/r8_00.prep.json"
for dt in DBNRec:r1_00 DBNRec:r2_00 DBNRec:r3_00 DBNRec:r4_00 DBNRec2:r5_00 DBNRec2:r6_00 DBNRec2:r7_00 DBNRec2:r4_09; do
  d=${dt%%:*}; t=${dt##*:}; pj=Riemann/Deep/$d/logs/$t.prep.json
  args=$(python3 -c "import json;from fractions import Fraction as F;p=json.loads(open('$pj').readline());print(p['row'],'%.9f'%float(F(p['u0'])),'%.9f'%float(F(p['u1'])),'%.9f'%float(F(p['T'])))")
  bash .scratch/memcap.sh 290 600 $UV Riemann/Deep/DBNRec2/code/prep_end.py $ARC $AFF $args $K/${d}_$t.dat > $K/${d}_$t.json 2> $K/${d}_$t.err
  if head -1 $K/${d}_$t.json | cmp -s - <(head -1 $pj); then log "regen $d/$t ($args) prep.json identical"; specs="$specs $K/${d}_$t.dat:$pj"; else log "regen $d/$t prep.json DIFFERS"; fi
done
step phase_arb 900 $UV $C/phase_arb.py $specs
rm -f $K/DBNRec*_r*.dat
step dense_r4 1200 $UV $C/falsify_dense.py 4
step scan_r8_32off 1500 $UV $C/scan_wrap.py 8 300 32
log "queue5 done"
