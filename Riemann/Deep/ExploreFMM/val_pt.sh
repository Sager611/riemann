#!/bin/bash
# POINT-mode (backward-error) certificates on the stored particle data: nu = m(1 - 1e-9) must pass, m(1 + 1e-9) must fail
cd "$(dirname "$0")"
for spec in "$@"; do
  set -- $spec; f=$1; P=$2; m=$3
  lo=$(python3 -c "print('%.12e' % ($m*(1-1e-9)))"); hi=$(python3 -c "print('%.12e' % ($m*(1+1e-9)))")
  echo "== $f $P margin $m"
  SLT_POINT=1 /usr/bin/time -l ./slt_pt runs/$f $P -l runs/$f.gen $lo $hi 2> runs/$f.pt.time | tee runs/$f.pt.json | python3 -c "
import sys, json
for l in sys.stdin:
    d = json.loads(l); print(d['nu'], 'ok' if d['inertia_ok'] else 'FAIL', 'U', d['U'], 'B', d['B'], 'kappa_eta %.2e' % d['kappa_eta'], 'eta %.2e' % d['eta'], 'kappa %.2e' % d['kappa'], 't', d['time'])"
  rg -o '[0-9]+  maximum resident set size' runs/$f.pt.time
done
