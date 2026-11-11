# self-tests: W two ways, W by mpmath contour quadrature, M1 d-sum vs direct double sum
import numpy as np, mpmath
from mollifier import Wfun, coeffs, main_terms
from math import gcd

for sigma in [0.52, 0.6, 0.75]:
    W = Wfun(sigma)
    for xi in [1.0, 1.5, 3.0, 10.0]:
        a = W.series(np.array([xi]))[0]
        b = W.asym(np.array([xi]))[0]
        mpmath.mp.dps = 20
        f = lambda u: (mpmath.gamma(1 + (0.5 + 1j * u)) * mpmath.zeta(2 * sigma + 0.5 + 1j * u)
                       / mpmath.zeta(2 * sigma) * mpmath.zeta(1.5 + 1j * u) * mpmath.mpf(xi) ** (0.5 + 1j * u))
        c = mpmath.quad(f, [-60, -20, 0, 20, 60]) / (2 * mpmath.pi)
        print("sigma=%.2f xi=%5.1f series=%.15f asym=%.15f contour=%.12f" % (sigma, xi, a, b, float(mpmath.re(c))))

# direct double-sum check of M1 and S2
y, sigma = 12.0, 0.6
N = int(y ** 1.5)
a, mu, primes, W = coeffs(y, sigma, N)
M1, S2 = main_terms(a, mu, primes, sigma)
d1 = d2 = 0.0
for h in range(1, N + 1):
    for k in range(1, N + 1):
        if a[h] and a[k]:
            g = gcd(h, k)
            l = h * k // g
            d1 += a[h] * a[k] / l ** (2 * sigma)
            d2 += a[h] * a[k] * g ** (2 - 2 * sigma) / (h * k)
d1 *= float(mpmath.zeta(2 * sigma))
print("M1 dsum=%.15f direct=%.15f | S2 dsum=%.15f direct=%.15f" % (M1, d1, S2, d2))
