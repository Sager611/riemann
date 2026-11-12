"""Float landscape of the complex (= twisted) pencil: sup_v lambda_max/kappa^2 for r in list, kappa in list; plus closure
ceilings for r = 4, 5. usage: cland.py r_list kappa_list N [ceil_r_list]"""
import sys, time, resource, json
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, HERE)
from ctest import mats, best
t0 = time.time()
rs = [int(x) for x in sys.argv[1].split(',')]; ks = [float(x) for x in sys.argv[2].split(',')]; N = int(sys.argv[3])
cr = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else []
out = {}
for r in rs:
    for kap in ks:
        D = 36 if kap < 5 else 44
        res = []
        for cp in (True, False):
            M0, MS, MSS, _ = mats(kap * np.pi, 0.25, r, N, D, cp)
            lam, v = best(M0, MS, MSS)
            res.append((lam / kap ** 2, v))
        out['r%d_k%s' % (r, kap)] = res
        print('r=%d kappa=%.3f: complex %.7f (v=%.4f) | real %.7f (v=%.4f) [%.0fs]' % (r, kap, res[0][0], res[0][1], res[1][0], res[1][1], time.time() - t0), flush=True)
for r in cr:
    f = lambda k: best(*mats(k * np.pi, 0.25, r, N, 36, True)[:3])[0] / k ** 2 - 1
    k = brentq(f, 2.8, 3.4, xtol=2e-6)
    fr = lambda k: best(*mats(k * np.pi, 0.25, r, N, 36, False)[:3])[0] / k ** 2 - 1
    kr = brentq(fr, 2.8, 3.4, xtol=2e-6)
    out['ceil_r%d' % r] = (k, kr)
    print('r=%d N=%d: complex closure ceiling %.6f | real %.6f [%.0fs]' % (r, N, k, kr, time.time() - t0), flush=True)
json.dump(out, open(HERE + '/logs/cland.json', 'w'), indent=1)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
