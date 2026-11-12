# Check the numerical inputs of Chourasiya-Simonic's Cor 2.3 proof (TeX l. 371): 7.031e5 < M2(3000) < 7.032e5, 1.249e6 < M2(3000,6000) < 1.25e6,
# M2(A,B) = int_A^B |zeta(1/2+it)|^4 dt. zeta from python-flint (Arb, 64-bit balls); Gauss-Legendre (numpy nodes) on panels of length 0.5;
# two node counts (6, 10) give the quadrature error estimate (not a rigorous enclosure).
import numpy as np, sys, time
from flint import acb, ctx
ctx.prec = 64
def M2(A, B, nn, h=0.5):
    x, w = np.polynomial.legendre.leggauss(nn)
    tot = 0.0; a = float(A)
    while a < B - 1e-12:
        b = min(a + h, float(B)); c, r = (a + b)/2, (b - a)/2
        s = 0.0
        for xi, wi in zip(x, w):
            z = acb(0.5, c + r*xi).zeta(); m = float(abs(z).mid())
            s += wi*m**4
        tot += r*s; a = b
    return tot
t0 = time.time()
for nn in (6, 10):
    v1 = M2(0, 3000, nn); v2 = M2(3000, 6000, nn)
    print('GL %2d nodes/0.5: M2(3000) = %.4f   M2(3000,6000) = %.4f   [%.0f s]' % (nn, v1, v2, time.time() - t0)); sys.stdout.flush()
print('claims: 7.031e5 < M2(3000) < 7.032e5 ->', 7.031e5 < v1 < 7.032e5, '; 1.249e6 < M2(3000,6000) < 1.25e6 ->', 1.249e6 < v2 < 1.25e6)
