"""Float evidence for the open pieces (not a certificate).
(a) r = 2 closure ceiling vs Galerkin size N (corner modes would grow with N); (b) r = 2 landscape sup_v lam/kappa^2;
(c) r = 5, 6 closure ceilings. usage: evid.py a|b|c"""
import sys, time, resource
import numpy as np
from scipy.optimize import brentq
import gcore as G
import gcore_r as GR
th = 0.25; t0 = time.time(); part = sys.argv[1]
def sup_ratio(r, k, N, D=36):
    M0, MS, MSS, B = GR.mats_r(k * np.pi, th, r, N, D)
    lam, v = G.best_v(M0, MS, MSS, lo=0.6, hi=2.2)
    return lam / k ** 2, v
if part == 'a':
    for N in (8, 20, 40, 60):
        k = brentq(lambda kk: sup_ratio(2, kk, N)[0] - 1, 3.25, 3.35, xtol=1e-7)
        print('r=2 N=%d: closure ceiling %.7f [%.0fs]' % (N, k, time.time() - t0), flush=True)
if part == 'b':
    for k in (3.31, 3.35, 3.4, 3.6, 4.0, 5.0, 6.0, 8.0):
        s, v = sup_ratio(2, k, 16, 36 if k < 5 else 44)
        print('r=2 kappa=%.2f: sup_v lam/k^2 = %.6f (v=%.4f) [%.0fs]' % (k, s, v, time.time() - t0), flush=True)
if part == 'c':
    for r in (5, 6):
        for N in (8, 14):
            k = brentq(lambda kk: sup_ratio(r, kk, N)[0] - 1, 2.4, 3.4, xtol=1e-6)
            print('r=%d N=%d: closure ceiling %.6f [%.0fs]' % (r, N, k, time.time() - t0), flush=True)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
