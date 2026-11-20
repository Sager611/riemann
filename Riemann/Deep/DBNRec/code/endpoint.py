#!/usr/bin/env python3
# dbnhead endpoint.py: float chain endpoint (DBNHigh chain.py) with certified window credits on some rows and ceilaudit5's model credits elsewhere.
# credit_i = min(tri_i - Scert_i, cr*_i) on certified rows (far shells u > u*: model), cr*_i elsewhere; dS_i = credit_i - max(0, ray_i).
# usage: endpoint.py ARC selfcons.json certs.json   (certs: {"row": [Scert, tri], ...})
import sys, json
arc, scf, cf = sys.argv[1:4]
res = json.load(open(scf)); certs = {int(k): tuple(v) for k, v in json.load(open(cf)).items()}   # row: [Scert, tri]
sys.argv = [sys.argv[0], arc, '/dev/null']
src = open('Riemann/Deep/DBNHigh/chain.py').read(); src = src[:src.index('t0 = time.time()')]
exec(src)
E0 = run()[0]; print('base endpoint %r' % E0, flush=True)
for kap in ('1.32', '7'):
    cr = [r['cr_' + kap] for r in res]
    E1 = run(dS=[max(0.0, c - max(0.0, r['ray'])) for c, r in zip(cr, res)])[0]
    cr2 = list(cr)
    print('kappa', kap, 'model endpoint %r (gain %.3e)' % (E1, E0 - E1), flush=True)
    for i, (S, tri) in certs.items():
        cw = tri - S; cr2[i] = min(cw, cr[i]); print('  row %d certified window credit %.6f model cr* %.6f -> used %.6f' % (i, cw, cr[i], cr2[i]))
    E2 = run(dS=[max(0.0, c - max(0.0, r['ray'])) for c, r in zip(cr2, res)])[0]
    print('kappa', kap, 'certified rows + model elsewhere: endpoint %r (gain %.3e)' % (E2, E0 - E2), flush=True)
