#!/usr/bin/env python3
# dbnrec3 credits3.py (tables only, no computation): credit table from cert_r<row>.json (rows 1-3: tranche 1, rows 4-7: tranche 2, read-only; rows 8-19: DBNRec3/logs).
# usage: credits3.py row ...   -> markdown table on stdout, JSON (exact A', B') to Riemann/Deep/DBNRec3/logs/credits.json
import sys, json, os
out = []
print('| row | source | A′ | B′ | Raymond A | Raymond B | dA | dB | shells | u_F |'); print('|---|---|---|---|---|---|---|---|---|---|')
for r in map(int, sys.argv[1:]):
    f, src = (('Riemann/Deep/DBNRec/logs/cert_r%d.json', 'tranche 1') if r <= 3 else ('Riemann/Deep/DBNRec2/logs/cert_r%d.json', 'tranche 2') if r <= 7 else ('Riemann/Deep/DBNRec3/logs/cert_r%d.json', 'tranche 3'))
    f = f % r
    if not os.path.exists(f): continue
    c = json.load(open(f)); A0 = c['A_raymond']; B0 = c['B_raymond']
    def fl(q): n, d = (q.split('/') + ['1'])[:2]; return int(n)/int(d)
    print('| %d | %s | %.9f | %.9f | %.9f | %.9f | %.6f | %.6f | %d | %s |' % (r, src, c['A_new_float'], c['B_new_float'], fl(A0), fl(B0), c['dA'], c['dB'], len(c['shells']), c['far']['uF']))
    out.append(dict(row=r, source=src, A_new=c['A_new'], B_new=c['B_new'], A_raymond=A0, B_raymond=B0, dA=c['dA'], dB=c['dB'], shells=len(c['shells']), uF=c['far']['uF']))
json.dump(out, open('Riemann/Deep/DBNRec3/logs/credits.json', 'w'), indent=1)
