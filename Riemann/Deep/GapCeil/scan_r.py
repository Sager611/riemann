"""r = 1,2,3: L2 closure ceiling and poly-degree values (float)."""
import sys, time, resource, json
import numpy as np
from scipy.optimize import brentq
import gcore as G
import gcore_r as GR
th = 0.25
t0 = time.time()
rs = [int(x) for x in sys.argv[1].split(',')]
res = {}
for r in rs:
    for N in (8, 14, 20):
        cache = {}
        def ratio(k):
            if k not in cache:
                M0, MS, MSS, B = GR.mats_r(k * np.pi, th, r, N, 36)
                cache[k] = G.best_v(M0, MS, MSS, lo=0.6, hi=2.2)
            return cache[k][0] / k ** 2
        lo, hi = 2.6, 3.6
        if ratio(lo) < 1: lo = 2.0
        k = brentq(lambda kk: ratio(kk) - 1, lo, hi, xtol=1e-6)
        res['r%d_N%d' % (r, N)] = k
        print('r=%d N=%d: L2-closure ceiling kappa = %.6f  (v=%.4f) [%.0fs]' % (r, N, k, cache[min(cache, key=lambda x: abs(x - k))][1], time.time() - t0), flush=True)
json.dump(res, open('logs/scan_r.json', 'w'), indent=1)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
