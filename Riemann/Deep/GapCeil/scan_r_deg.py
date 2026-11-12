"""r = 2: kappa_max for polynomial P of degree d (phi = I^{2r} P_d subspace) vs the L2 closure (float)."""
import sys, time, resource, json
import numpy as np
from scipy.optimize import brentq
import gcore as G
import gcore_r as GR
th = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
r = int(sys.argv[1]); degs = [int(x) for x in sys.argv[3].split(',')]
t0 = time.time(); Nb = max(degs) + 2 * r + 4
cache = {}
def mats(k):
    if k not in cache:
        cache[k] = GR.mats_r(k * np.pi, th, r, Nb, 36)
    return cache[k]
res = {}
for d in degs:
    def ratio(k):
        M0, MS, MSS, B = mats(k)
        Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
        Q, _ = np.linalg.qr(Ir)
        m0, ms, mss = (Q.T @ A @ Q for A in (M0, MS, MSS))
        return G.best_v(m0, ms, mss, lo=0.6, hi=2.2)[0] / k ** 2
    try:
        k = brentq(lambda kk: ratio(kk) - 1, 2.8, 3.4, xtol=1e-6)
    except Exception as e:
        k = float('nan')
    res[d] = k
    print('r=%d theta=%.4f deg %2d: kappa_max = %.6f [%.0fs]' % (r, th, d, k, time.time() - t0), flush=True)
json.dump(res, open('logs/scan_r%d_deg_th%s.json' % (r, th), 'w'), indent=1)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
