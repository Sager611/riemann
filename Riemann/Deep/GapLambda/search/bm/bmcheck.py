import sys, time, resource
sys.path.insert(0, '.')
import numpy as np
from numpy.polynomial import Polynomial, Legendre
import bmfloat as bf
t0 = time.time()
PLEG8 = np.array([1., 2.8691, 4.6865, 5.6252, 7.1055, 6.0975, 6.9685, 3.1024, 5.6301])
def best_v(M, k, P=None):
    vs = np.linspace(1.0, 1.6, 61)
    f = (lambda v: bf.h1max(M, k, v)[0]) if P is None else (lambda v: bf.h1P(M, k, v, P))
    hs = [f(v) for v in vs]; i = int(np.argmax(hs)); v, h, dv = vs[i], hs[i], 0.005
    while dv > 1e-7:
        moved = False
        for vv in (v - dv, v + dv):
            hh = f(vv)
            if hh > h: h, v, moved = hh, vv, True
        if not moved: dv /= 2
    return h, v
# 1) quadrature convergence at kappa = 3.25 (deg 8): free-P optimum and the fixed n=8 optimum P
for n in (8, 10, 12, 14):
    M, im = bf.mats(0.25, 1, 3.25, 9, n)
    h, v = best_v(M, 3.25); hp, vp = best_v(M, 3.25, PLEG8)
    print('n=%2d kappa=3.25 deg8: max_P h1 = %.8f (v=%.5f) | fixed P: h1 = %.8f (v=%.5f) imag %.0e (%.0fs)' % (n, h, v, hp, vp, im, time.time() - t0), flush=True)
# 2) theta strictly below 1/4
for th in (0.2499, 0.249):
    M, im = bf.mats(th, 1, 3.25, 9, 10)
    h, v = best_v(M, 3.25)
    print('theta=%.4f n=10 kappa=3.25: max_P h1 = %.8f (v=%.5f)' % (th, h, v), flush=True)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
