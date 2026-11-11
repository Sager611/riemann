"""Amplifier exponent r = 2, 3: certified (Arb) bounds for P = 1 and the best r = 2 polynomial found."""
import sys
sys.path.insert(0, '.')
import numpy as np
from flint import fmpq
import cert


def best_uv(th, r, P):
    mt = cert.main_terms(th, r, P)
    c = float(mt['c'][0].mid()); p1 = [float(x.mid()) for x in mt['c1']]; p2 = [float(x.mid()) for x in mt['c2']]
    best = (-1, 0, 0)
    for u in np.linspace(0.3, 1.0, 701):
        C1 = np.polyval(p1[::-1], u); C2 = np.polyval(p2[::-1], u)
        vs = np.linspace(0, 1, 1001); l0 = (1 + 4 * vs + np.sqrt(1 + 8 * vs)) / 8
        k2 = (np.sqrt(144 * vs ** 2 * C2 ** 2 + 48 * l0 * c * C1) - 12 * vs * C2) / (16 * C1)
        i = int(np.argmax(k2))
        if k2[i] > best[0]: best = (k2[i], u, vs[i])
    return fmpq(round(best[1] * 10000), 10000), fmpq(round(best[2] * 10000), 10000)


th = fmpq(1249, 10000)
for r in (2, 3):
    u, v = best_uv(th, r, [fmpq(1)])
    cert.bound(th, r, [fmpq(1)], u, v, label='r=%d P=1' % r)
pl = ['1', '1.0376', '-1.4821', '-2.3091', '1.7174', '6.098', '6.7703']
P = [fmpq(0)] * len(pl)
for m, s in enumerate(pl):
    for k, cf in enumerate(cert.shleg_mono(m)):
        neg = s.startswith('-'); t = s.lstrip('-'); a, b = (t.split('.') + [''])[:2]
        q = fmpq(int(a + b), 10 ** len(b)); q = -q if neg else q
        P[k] += q * cf
u, v = best_uv(th, 2, P)
cert.bound(th, 2, P, u, v, label='r=2 best deg-6 P from search')
