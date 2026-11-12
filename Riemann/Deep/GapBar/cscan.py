"""Complex-amplifier (Hermitian pencil) L2-closure ceilings, float. usage: cscan.py r_list N"""
import sys, time, resource, json
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, HERE)
from ctest import mats, best
t0 = time.time()
rs = [int(x) for x in sys.argv[1].split(',')]; N = int(sys.argv[2])
res = {}
for r in rs:
    cache = {}
    def ratio(k):
        if k not in cache:
            M0, MS, MSS, _ = mats(k * np.pi, 0.25, r, N, 36, True)
            cache[k] = best(M0, MS, MSS)
            print('  r=%d k=%.6f ratio=%.7f v=%.4f [%.0fs]' % (r, k, cache[k][0] / k ** 2, cache[k][1], time.time() - t0), flush=True)
        return cache[k][0] / k ** 2
    k = brentq(lambda kk: ratio(kk) - 1, 3.1, 3.6, xtol=2e-6)
    res['r%d_N%d' % (r, N)] = k
    print('r=%d N=%d: complex closure ceiling kappa = %.6f [%.0fs]' % (r, N, k, time.time() - t0), flush=True)
json.dump(res, open(HERE + '/logs/cscan_%s_N%d.json' % (sys.argv[1].replace(',', ''), N), 'w'), indent=1)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
