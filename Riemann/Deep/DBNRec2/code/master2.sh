#!/bin/bash
# dbnrec2 master2.sh: one sequential queue (at most one compute process at a time; every step under memcap 290 MB), functions as tranche 1's master.sh.
# Row 4: shells 06-10 at tranche 1's T = 0.834173 and shell list (shells 00-05: tranche 1's logs, copied; .dat regenerated, prep lines identical), certify 00-10.
# Then rows 5, 6, ... (dorow: shell 0 at T0, T from its dump, outer shells), each certified. Before each new row: optional hook (.scratch/dbnrec2/hook.sh).
# NEWROW_DEADLINE / SHELL_DEADLINE (epoch s) stop new rows / new shells (re-read from .scratch/dbnrec2/deadlines before each row; rows with a cert are skipped). The exact chain is run separately (chain2.sh).
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=Riemann/Deep/DBNRec2/logs; C=Riemann/Deep/DBNRec2/code; W=.scratch/dbnrec2/ws
ARC=.scratch/dbnrec2/dl/Proof_and_Replay__Lambda_Peer_Review__verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json
AFF=.scratch/dbnrec2/dl/Proof_and_Replay__Lambda_Peer_Review__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
UVF="uv run --no-project --with mpmath --with python-flint"
NEWROW_DEADLINE=${NEWROW_DEADLINE:-0}; export SHELL_DEADLINE
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/master.log; }
cert() { r=$1; tags=$(ls $L/r${r}_*.bb.out 2>/dev/null | sed 's#.*/##; s#.bb.out##' | sort | tr '\n' ' ')
  /usr/bin/time -l bash .scratch/memcap.sh 290 1200 $UVF python3 $C/certify.py $r $tags > $L/cert_r$r.json 2> $L/cert_r$r.err
  log "cert row $r rc=$? $(python3 -c "import json;c=json.load(open('$L/cert_r$r.json'));print('A',c['A_new_float'],'B',c['B_new_float'],'dA',c['dA'],'dB',c['dB'])" 2>&1 | tail -1)"; }
dorow() { r=$1
  bash .scratch/memcap.sh 290 600 $UVF python3 $C/prep_end.py $ARC $AFF $r 0 0.01 0.7 .scratch/dbnrec2/tmp.dat > .scratch/dbnrec2/tmp.json 2>/dev/null
  T0=$(python3 $C/plan_row.py t0 $r .scratch/dbnrec2/tmp.json); log "row $r T0 $T0"
  bash $C/run_row.sh $r $T0 0 0.01
  if ! rg -q '^CERT' $L/r${r}_00.bb.out; then T0=$(python3 -c "print('%.6f'%($T0+0.003))"); rm -f $L/r${r}_00.*; log "row $r retry T0 $T0"; bash $C/run_row.sh $r $T0 0 0.01; fi
  T=$(python3 $C/plan_row.py T $r $L/r${r}_00.prep.json $L/r${r}_00.unres.txt)
  SH=$(python3 $C/plan_row.py shells $L/r${r}_00.prep.json $T); log "row $r T $T shells $SH"
  bash $C/run_row.sh $r $T 0 $SH
  cert $r; }
hook() { if [ -f .scratch/dbnrec2/hook.sh ]; then mv .scratch/dbnrec2/hook.sh .scratch/dbnrec2/hook.run; log "hook start"; bash .scratch/dbnrec2/hook.run >> $L/hook.log 2>&1; log "hook done rc=$?"; fi; }
log "master2 start (NEWROW_DEADLINE $NEWROW_DEADLINE SHELL_DEADLINE $SHELL_DEADLINE)"
$UVF python3 $C/check.py 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 > $L/check.out 2>&1; log "check.py rc=$? $(tail -1 $L/check.out | cut -c1-200)"
log "row 4 T 0.834173 shells (tranche 1) 0.01 0.02 0.035 0.05 0.07 0.09 0.11 0.13 0.15 0.17 0.1714"
if [ ! -s $L/cert_r4.json ]; then bash $C/run_row.sh 4 0.834173 0 0.01 0.02 0.035 0.05 0.07 0.09 0.11 0.13 0.15 0.17 0.1714; cert 4; fi
for r in 5 6 7 8 9 10 11 12 13 14 15; do
  hook
  [ -f .scratch/dbnrec2/deadlines ] && read NEWROW_DEADLINE SHELL_DEADLINE < .scratch/dbnrec2/deadlines && export SHELL_DEADLINE   # re-read before each row
  [ -s $L/cert_r$r.json ] && continue
  [ $(date +%s) -gt $NEWROW_DEADLINE ] && { log "new-row deadline: stop before row $r"; break; }
  [ $(date +%s) -gt $SHELL_DEADLINE ] && { log "shell deadline: stop before row $r"; break; }
  dorow $r
done
log "master2 done"
