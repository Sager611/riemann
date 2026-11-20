#!/bin/bash
# dbnrecaudit queue9.sh (third reviewer pass, 17:43 UTC: "B7: phase-error bounds for 77 unchecked shells"): after queue7.sh (edges), regenerate
# the .dat of EVERY logged shell of rows 1-8 (86; prep.json line must equal the logged one) and run phase_arb.py on all of them.
# One process at a time, memcap 290 MB; .dat files in /tmp/dbnra_all (≈ 70 MB), deleted at the end.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=Riemann/Deep/DBNRecAudit/logs; C=Riemann/Deep/DBNRecAudit/code; K=/tmp/dbnra_all; mkdir -p $K
UV="uv run --no-project --with mpmath --with numpy --with python-flint python3"
DL=.scratch/dbnrecaudit/dl/Proof_and_Replay__Lambda_Peer_Review__verification
ARC=${DL}__arcs__ACTUAL_67_ROW_ARC_INPUTS.json; AFF=${DL}__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/queue.log; }
while pgrep -f 'DBNRecAudit/code/queue7.sh' > /dev/null; do sleep 3; done
log "start regen_all"
specs=""; nid=0; ndiff=0
for dt in $(python3 -c "import json;print(' '.join(s['d']+':'+s['tag'] for s in json.load(open('$L/check_logs_shells.json'))))"); do
  d=${dt%%:*}; t=${dt##*:}; pj=Riemann/Deep/$d/logs/$t.prep.json
  args=$(python3 -c "import json;from fractions import Fraction as F;p=json.loads(open('$pj').readline());print(p['row'],'%.9f'%float(F(p['u0'])),'%.9f'%float(F(p['u1'])),'%.9f'%float(F(p['T'])))")
  bash .scratch/memcap.sh 290 600 $UV Riemann/Deep/DBNRec2/code/prep_end.py $ARC $AFF $args $K/${d}_$t.dat > $K/${d}_$t.json 2> /dev/null
  if head -1 $K/${d}_$t.json | cmp -s - <(head -1 $pj); then nid=$((nid+1)); specs="$specs $K/${d}_$t.dat:$pj"; else ndiff=$((ndiff+1)); log "regen_all $d/$t prep.json DIFFERS"; fi
done
log "end regen_all: $nid identical, $ndiff differ"
/usr/bin/time -l bash .scratch/memcap.sh 290 900 $UV $C/phase_arb.py $specs > $L/phase_arb_all.out 2> $L/phase_arb_all.err
log "end phase_arb_all rc=$? $(rg -o 'memcap: peak [0-9]+ MB.*' $L/phase_arb_all.err | tail -1) $(rg -o '[0-9]+  maximum resident' $L/phase_arb_all.err | tail -1)"
rm -rf $K
log "queue9 done"
