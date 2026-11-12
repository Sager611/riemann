"""ceiling scan: sup_v lambda_max(kappa)/kappa^2 for the L2(phi) pencil (closure) and for deg-d polynomial P."""
import sys, time, resource, json
import numpy as np
from scipy.optimize import brentq
import gcore as G

th = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
Ns = [int(s) for s in (sys.argv[2] if len(sys.argv) > 2 else '8,12,16,24,32,40').split(',')]
degs = [int(s) for s in (sys.argv[3] if len(sys.argv) > 3 else '2,8,14').split(',')]
D = 36
t0 = time.time()
cache = {}

def get(kap, N):
    key = round(kap, 12)
    if key not in cache:
        cache[key] = G.mats(kap * np.pi, th, max(Ns + [max(degs) + 4]), D)
    M0, MS, MSS, B = cache[key]
    return M0, MS, MSS, B

def ratio(kap, N=None, deg=None):
    M0, MS, MSS, B = get(kap, N)
    if deg is not None:
        Nb = M0.shape[0]
        I2 = (B.I @ B.I)[:Nb, :deg + 1]
        Q, _ = np.linalg.qr(I2)
        M0, MS, MSS = (Q.T @ A @ Q for A in (M0, MS, MSS))
    else:
        M0, MS, MSS = (A[:N, :N] for A in (M0, MS, MSS))
    lam, v = G.best_v(M0, MS, MSS)
    return lam / kap ** 2, v

res = {}
for deg in degs:
    f = lambda k: ratio(k, deg=deg)[0] - 1
    k = brentq(f, 3.0, 3.6, xtol=1e-7)
    res['deg%d' % deg] = (k, ratio(k, deg=deg)[1])
    print('poly deg %2d: kappa_max = %.6f (v=%.5f)  [%.0fs]' % (deg, k, res['deg%d' % deg][1], time.time() - t0), flush=True)
for N in Ns:
    f = lambda k: ratio(k, N=N)[0] - 1
    k = brentq(f, 3.0, 4.5, xtol=1e-7)
    res['N%d' % N] = (k, ratio(k, N=N)[1])
    print('L2(phi) N=%2d: kappa_max = %.6f (v=%.5f)  [%.0fs]' % (N, k, res['N%d' % N][1], time.time() - t0), flush=True)
json.dump(res, open('scan_th%s.json' % th, 'w'), indent=1)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
