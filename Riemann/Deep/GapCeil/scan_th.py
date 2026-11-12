"""L2-closure ceiling as a function of the mollifier length theta (same main-term formula; BBLR proves it only for theta < 1/4)."""
import sys, time, resource, json
import numpy as np
from scipy.optimize import brentq
import gcore as G
t0 = time.time(); res = {}
for th in [float(x) for x in sys.argv[1].split(',')]:
    cache = {}
    def ratio(k):
        if k not in cache:
            M0, MS, MSS, B = G.mats(k * np.pi, th, 16, 40)
            cache[k] = G.best_v(M0, MS, MSS, lo=0.6, hi=2.4)
        return cache[k][0] / k ** 2
    k = brentq(lambda kk: ratio(kk) - 1, 2.0, 4.5, xtol=1e-6)
    res[th] = k
    print('theta=%.4f: L2-closure ceiling kappa* = %.6f [%.0fs]' % (th, k, time.time() - t0), flush=True)
json.dump(res, open('logs/scan_th.json', 'w'), indent=1)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
