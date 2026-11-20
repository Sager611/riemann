#!/bin/bash
# dbnrec master.sh: one sequential queue (at most one compute process at a time; every step under memcap 290 MB).
# After row 1's shells: certify row 1, exact chain (producer + verifier) with row 1; then rows 2..7 (shell 0 at T0, T from its dump, outer shells),
# each certified; final exact chain with all certified rows. NEWROW_DEADLINE / SHELL_DEADLINE (epoch s) stop new work in time for the write-up.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=Riemann/Deep/DBNRec2/logs; C=Riemann/Deep/DBNRec2/code; W=.scratch/dbnrec2/ws
ARC=.scratch/dbnrec2/dl/Proof_and_Replay__Lambda_Peer_Review__verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json
AFF=.scratch/dbnrec2/dl/Proof_and_Replay__Lambda_Peer_Review__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
UVF="uv run --no-project --with mpmath --with python-flint"
UVP="uv run --no-project --python 3.12 --with python-flint==0.8.0"
NEWROW_DEADLINE=${NEWROW_DEADLINE:-0}
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/master.log; }
cert() { r=$1; tags=$(ls $L/r${r}_*.bb.out 2>/dev/null | sed 's#.*/##; s#.bb.out##' | sort | tr '\n' ' ')
  /usr/bin/time -l bash .scratch/memcap.sh 290 1200 $UVF python3 $C/certify.py $r $tags > $L/cert_r$r.json 2> $L/cert_r$r.err
  log "cert row $r rc=$? $(python3 -c "import json;c=json.load(open('$L/cert_r$r.json'));print('A',c['A_new_float'],'B',c['B_new_float'],'dA',c['dA'],'dB',c['dB'])" 2>&1 | tail -1)"; }
chain() { tag=$1; shift; certs=""; for r in "$@"; do certs="$certs $L/cert_r$r.json"; done
  $UVF python3 $C/mkrec.py source $certs > $L/mkrec_$tag.json 2>&1; log "mkrec source $tag: $(cut -c1-300 $L/mkrec_$tag.json)"
  rm -f $W/dbnrec/FLOW_REC.json.gz
  /usr/bin/time -l bash .scratch/memcap.sh 290 1800 $UVP python3 $C/make_flow_rec.py --workspace $W --out $W/dbnrec/FLOW_REC.json.gz > $L/make_flow_$tag.log 2>&1
  log "make_flow $tag rc=$? $(rg -o '\"endpoint_decimal\": [0-9.e-]+' $L/make_flow_$tag.log)"
  cp $W/dbnrec/FLOW_REC_SUMMARY.json $L/FLOW_REC_SUMMARY_$tag.json
  $UVF python3 $C/mkrec.py verifier > $L/mkrec_ver_$tag.json 2>&1; log "mkrec verifier $tag: $(cut -c1-200 $L/mkrec_ver_$tag.json)"
  (cd .scratch/dbnrec2/verify && /usr/bin/time -l bash ../../memcap.sh 290 1800 $UVP python3 verify_arcs_rec.py --workspace ../ws > ../../../$L/verify_$tag.log 2>&1)
  log "verify $tag rc=$? $(rg -o '\"PASS\": [a-z]+' $L/verify_$tag.log)"
  cp .scratch/dbnrec2/verify/EXACT_ARC_VERIFICATION_FRESH_V1.json $L/EXACT_ARC_VERIFICATION_$tag.json 2>/dev/null
  rm -f .scratch/dbnrec2/verify/ALL_221010_ARC_INEQUALITY_WITNESSES.jsonl.gz; }
dorow() { r=$1
  bash .scratch/memcap.sh 290 600 $UVF python3 $C/prep_end.py $ARC $AFF $r 0 0.01 0.7 .scratch/dbnrec2/tmp.dat > .scratch/dbnrec2/tmp.json 2>/dev/null
  T0=$(python3 $C/plan_row.py t0 $r .scratch/dbnrec2/tmp.json); log "row $r T0 $T0"
  bash $C/run_row.sh $r $T0 0 0.01
  if ! rg -q '^CERT' $L/r${r}_00.bb.out; then T0=$(python3 -c "print('%.6f'%($T0+0.003))"); rm -f $L/r${r}_00.*; log "row $r retry T0 $T0"; bash $C/run_row.sh $r $T0 0 0.01; fi
  T=$(python3 $C/plan_row.py T $r $L/r${r}_00.prep.json $L/r${r}_00.unres.txt)
  SH=$(python3 $C/plan_row.py shells $L/r${r}_00.prep.json $T); log "row $r T $T shells $SH"
  bash $C/run_row.sh $r $T 0 $SH
  cert $r; }
while pgrep -f 'run_row.sh 1 ' > /dev/null; do sleep 10; done
log "row 1 shells done"
cert 1
chain A 1
done_rows="1"
for r in 2 3 4 5 6 7; do
  [ $(date +%s) -gt $NEWROW_DEADLINE ] && { log "new-row deadline: stop before row $r"; break; }
  dorow $r; done_rows="$done_rows $r"
done
chain B $done_rows
log "master done rows: $done_rows"
