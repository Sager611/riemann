import sys, time
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, HERE)
from ctest import mats, best
t0 = time.time()
for r in (4, 5):
    for N in (8, 10):
        for cp in (True, False):
            try:
                f = lambda k: best(*mats(k * np.pi, 0.25, r, N, 36, cp)[:3])[0] / k ** 2 - 1
                k = brentq(f, 2.7, 3.4, xtol=2e-6)
                print('r=%d N=%d %s closure ceiling %.6f [%.0fs]' % (r, N, 'complex' if cp else 'real', k, time.time() - t0), flush=True)
            except Exception as e:
                print('r=%d N=%d %s failed: %s' % (r, N, 'complex' if cp else 'real', str(e)[:80]), flush=True)
