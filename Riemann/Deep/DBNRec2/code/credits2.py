#!/usr/bin/env python3
# dbnrec2 credits2.py (tables only, no computation): credit table from cert_r<row>.json (rows 1-3: tranche 1 read-only; others: DBNRec2/logs).
# usage: credits2.py row ...   -> markdown table on stdout, JSON (exact A', B') to Riemann/Deep/DBNRec2/logs/credits.json
import sys, json
from fractions import Fraction as Fr
out = []; print('| row | source | A′ | B′ | Raymond A | Raymond B | dA | dB | shells | u_F |'); print('|---|---|---|---|---|---|---|---|---|---|')
for r in map(int, sys.argv[1:]):
    f = ('Riemann/Deep/DBNRec/logs/cert_r%d.json' if r <= 3 else 'Riemann/Deep/DBNRec2/logs/cert_r%d.json') % r
    c = json.load(open(f)); A1, B1, A0, B0 = (Fr(c[k]) for k in ('A_new', 'B_new', 'A_raymond', 'B_raymond'))
    print('| %d | %s | %.9f | %.9f | %.9f | %.9f | %.6f | %.6f | %d | %s |' % (r, 'tranche 1' if r <= 3 else 'tranche 2', A1, B1, A0, B0, A0 - A1, B0 - B1, len(c['shells']), c['far']['uF']))
    out.append(dict(row=r, cert=f, A_new=c['A_new'], B_new=c['B_new'], A_raymond=c['A_raymond'], B_raymond=c['B_raymond'], dA=float(A0 - A1), dB=float(B0 - B1)))
json.dump(out, open('Riemann/Deep/DBNRec2/logs/credits.json', 'w'), indent=1)
