import sys, time, resource
sys.path.insert(0, '.')
import numpy as np
from numpy.polynomial import Polynomial, Legendre
import bmfloat as bf
t0 = time.time()
# BM check: theta=1/4, r=1, v=1.26, P = 1 - 5.8x + 6.4x^2, kappa = 3.18 -> h1 > 1.0002 (l.333-337)
Pm = Polynomial([1, -5.8, 6.4])
Pl = Pm.convert(kind=Legendre, domain=[0, 1]).coef
for n in (8, 10, 12):
    M, im = bf.mats(0.25, 1, 3.18, 3, n)
    print('n=%d  h1(BM P, v=1.26, k=3.18) = %.6f   rel imag %.1e  (%.0fs)' % (n, bf.h1P(M, 3.18, 1.26, Pl), im, time.time() - t0), flush=True)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
