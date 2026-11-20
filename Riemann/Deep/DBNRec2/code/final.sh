#!/bin/bash
# dbnrec2 final.sh: after master2.sh ends (or at the stop time STOP, epoch s: kill it, drop the unfinished shell, certify that row from its finished shells),
# run the exact chain (chain2.sh F) with every certified row: tranche 1's cert_r1-3.json + tranche 2's cert_r4.. .
cd /Users/adriansager/playground/riemann
L=Riemann/Deep/DBNRec2/logs; C=Riemann/Deep/DBNRec2/code; STOP=$1
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/master.log; }
tree() { local p; for p in "$@"; do echo $p; tree $(pgrep -P $p); done; }
while pgrep -f 'bash Riemann/Deep/DBNRec2/code/master2.sh' > /dev/null && [ $(date +%s) -lt $STOP ]; do sleep 5; done
M=$(pgrep -f 'bash Riemann/Deep/DBNRec2/code/master2.sh')
if [ -n "$M" ]; then
  kill -9 $(tree $M) 2>/dev/null; sleep 1; log "final: stop time reached, master2 killed"
  for f in $L/r*_*.bb.out; do rg -q '^CERT' $f || { t=${f%.bb.out}; log "final: drop unfinished shell $(basename $t)"; rm -f $t.bb.out $t.prep.json $t.unres.txt; }; done
  r=$(ls $L/r*_*.bb.out | sed 's#.*/r##; s#_.*##' | sort -n | tail -1)
  python3 -c "import json;json.load(open('$L/cert_r$r.json'))['A_new']" 2>/dev/null || { tags=$(ls $L/r${r}_*.bb.out | sed 's#.*/##; s#.bb.out##' | sort | tr '\n' ' ')
    /usr/bin/time -l bash .scratch/memcap.sh 290 1200 uv run --no-project --with mpmath --with python-flint python3 $C/certify.py $r $tags > $L/cert_r$r.json 2> $L/cert_r$r.err; log "final: cert row $r rc=$? (shells $tags)"; }
fi
certs="Riemann/Deep/DBNRec/logs/cert_r1.json Riemann/Deep/DBNRec/logs/cert_r2.json Riemann/Deep/DBNRec/logs/cert_r3.json"
for r in $(seq 4 30); do f=$L/cert_r$r.json; [ -s $f ] && python3 -c "import json;json.load(open('$f'))['A_new']" 2>/dev/null && certs="$certs $f"; done
log "final: chain F with $certs"
bash $C/chain2.sh F $certs
log "final: done"
