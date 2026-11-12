"""gapbar: float optimiser of a REAL amplifier P for the MIDPOINT mollifier M(s + i kappa pi/(2L)) (BM Lemma main with
alpha3 = -beta3 = i kappa pi/(2L): the kernel gets the factor y^{-alpha3(x1+x2)-beta3(x3+x4)} = e^{ic(r2-r1)/2}, i.e. the
twisted real kernel 4 sin(c(X1-X2+sigma)/2) sin(c(X3-X4+sigma)/2)).  Wraps GapCeil r2opt.py with the twisted rho_parts.
usage: optmid.py r theta d kappa [scan]"""
import sys, json
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
sys.path.insert(0, HERE)
import gcore as G
import gcore_r as GR
orig = G.rho_parts
def rho_tw(X1, X2, X3, X4, c):
    out = orig(X1, X2, X3, X4, c)
    tw = np.exp(1j * c * (X3 + X4 - X1 - X2) / 2)
    return [o * tw for o in out]
G.rho_parts = rho_tw
r = int(sys.argv[1]); th = float(sys.argv[2]); d = int(sys.argv[3]); kap = float(sys.argv[4])
def setup(k):
    Nb = d + 2 * r + 6
    M0, MS, MSS, B = GR.mats_r(k * np.pi, th, r, Nb, 40)
    Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
    Q, R = np.linalg.qr(Ir)
    return [Q.T @ A @ Q for A in (M0, MS, MSS)], R, B
if len(sys.argv) > 5 and sys.argv[5] == 'scan':
    from scipy.optimize import brentq
    f = lambda k: G.best_v(*setup(k)[0], lo=0.8, hi=2.0)[0] / k ** 2 - 1
    print('r=%d theta=%.4f deg=%d: midpoint real-P kappa_max = %.6f' % (r, th, d, brentq(f, 3.0, 3.42, xtol=1e-6)), flush=True)
    raise SystemExit
(m0, ms, mss), R, B = setup(kap)
lam, v = G.best_v(m0, ms, mss, lo=0.8, hi=2.0)
ev, y = G.top(m0, ms, mss, v)
a = np.linalg.solve(R, y)
An = a * B.norm[:d + 1]
An = An / An[0]
print('r=%d theta=%.4f deg=%d kappa=%.5f: float midpoint h1 = %.8f at v = %.6f' % (r, th, d, kap, lam / kap ** 2, v))
json.dump({'r': r, 'theta': sys.argv[2], 'd': d, 'kappa': sys.argv[4], 'v': '%.6f' % v, 'h1_float': lam / kap ** 2,
           'A': ['%.15e' % x for x in An], 'midpoint': True},
          open(HERE + '/logs/optmid_r%d_d%d_k%s.json' % (r, d, sys.argv[4]), 'w'), indent=1)
