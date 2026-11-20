#!/bin/bash
# dbnrec3 final3.sh STOP: after master3.sh ends (or at the stop time STOP, epoch s: kill it, drop the unfinished shell, certify that row from its finished shells),
# run the exact chain (chain3.sh G) with every certified row: tranche 1's cert_r1-3.json, tranche 2's cert_r4-7.json, tranche 3's cert_r8.. (row 8: tranche 2's
# partial cert only if tranche 3's is missing). Rows whose certified fields are not below Raymond's are skipped by mkrec.py.
cd /Users/adriansager/playground/riemann
L=Riemann/Deep/DBNRec3/logs; C=Riemann/Deep/DBNRec3/code; STOP=$1
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/master.log; }
tree() { local p; for p in "$@"; do echo $p; tree $(pgrep -P $p); done; }
ok() { python3 -c "import json;json.load(open('$1'))['A_new']" 2>/dev/null; }
while pgrep -f 'bash Riemann/Deep/DBNRec3/code/master3.sh' > /dev/null && [ $(date +%s) -lt $STOP ]; do sleep 5; done
M=$(pgrep -f 'bash Riemann/Deep/DBNRec3/code/master3.sh')
if [ -n "$M" ]; then
  kill -9 $(tree $M) 2>/dev/null; sleep 1; log "final: stop time reached, master3 killed"
  for f in $L/r*_*.bb.out; do rg -q '^CERT' $f || { t=${f%.bb.out}; log "final: drop unfinished shell $(basename $t)"; rm -f $t.bb.out $t.prep.json $t.unres.txt; }; done
  r=$(ls $L/r*_*.bb.out | sed 's#.*/r##; s#_.*##' | sort -n | tail -1)
  ok $L/cert_r$r.json || { tags=$(ls $L/r${r}_*.bb.out | sed 's#.*/##; s#.bb.out##' | sort | tr '\n' ' ')
    /usr/bin/time -l bash .scratch/memcap.sh 290 1200 uv run --no-project --with mpmath --with python-flint python3 $C/certify.py $r $tags > $L/cert_r$r.json 2> $L/cert_r$r.err; log "final: cert row $r rc=$? (shells $tags)"; }
fi
certs="Riemann/Deep/DBNRec/logs/cert_r1.json Riemann/Deep/DBNRec/logs/cert_r2.json Riemann/Deep/DBNRec/logs/cert_r3.json"
for r in 4 5 6 7; do certs="$certs Riemann/Deep/DBNRec2/logs/cert_r$r.json"; done
if ok $L/cert_r8.json; then certs="$certs $L/cert_r8.json"; else certs="$certs Riemann/Deep/DBNRec2/logs/cert_r8.json"; fi
for r in 9 10 11 12 13 14 15 17 18 19; do f=$L/cert_r$r.json; [ -s $f ] && ok $f && certs="$certs $f"; done
log "final: chain G with $certs"
bash $C/chain3.sh G $certs
log "final: done"
