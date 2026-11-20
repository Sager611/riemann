#!/usr/bin/env python3
# dbnrec endpoint_rec.py: float estimate (DBNHigh chain.py, scheme 1) of the chain endpoint with certified (A', B') on some rows.
# row source a + b (t - t-) with a = M - N - A, b = (A - B)/(t+ - t-)  ->  a' = a + (A - A'), b' = (A' - B')/(t+ - t-).
# usage: endpoint_rec.py ARC cert_r<i>.json ...
import sys, json
from fractions import Fraction as F
arc, certs = sys.argv[1], sys.argv[2:]
sys.argv = [sys.argv[0], arc, '/dev/null']
src = open('Riemann/Deep/DBNHigh/chain.py').read(); src = src[:src.index('t0 = time.time()')]
exec(src)
E0 = run()[0]; print('base endpoint %r' % E0, flush=True)
for cf in certs:
    c = json.load(open(cf)); i = c['row']; r = R[i]
    A0, B0, A1, B1 = (float(F(c[k])) for k in ('A_raymond', 'B_raymond', 'A_new', 'B_new'))
    if A1 <= A0 and B1 <= B0:
        r['a'] += A0 - A1; r['b'] = (A1 - B1)/(r['end'] - r['start'])
        print('row %d: dA %.6f dB %.6f' % (i, A0 - A1, B0 - B1), flush=True)
    E1 = run()[0]; print('  endpoint with rows so far %r (gain %.4e)' % (E1, E0 - E1), flush=True)
