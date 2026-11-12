# t9_consts.py (explore39): small constants for the report.
import math
from scipy.special import digamma
from scipy.optimize import brentq
D = float(digamma(0.25)) - math.log(math.pi)
K = lambda r: math.atan(math.exp(-r/2)) + 0.5*math.log(1/math.tanh(r/4))
kap = lambda z: math.exp(z/2)/(math.exp(z) - math.exp(-z))
J = lambda z: kap(z) - 2*math.cosh(z/2)
rho = brentq(lambda v: v**3 - v - 1, 1, 2)
print('D = %.10f, rho = %.12f, rho^2 = %.12f, ln rho = %.6f' % (D, rho, rho**2, math.log(rho)))
print('J(ln rho) = %.2e' % J(math.log(rho)))
# killing bound W(0) = D + 2K(a) + 8 sinh(a/2) (explore32 S2 potential at x = 0) -> zero
W0 = lambda a: D + 2*K(a) + 8*math.sinh(a/2)
a0 = brentq(W0, 0.05, 0.3); print('inf W = W(0) = 0 at a = %.6f, c = %.6f' % (a0, math.exp(2*a0)))
print('W(0) at c = 1.2: %.4f' % W0(math.log(1.2)/2))
# convexity check of J on (0, ln rho^2]
import numpy as np
zs = np.linspace(1e-3, 2*math.log(rho), 2000); h = 1e-4
J2 = [(J(z + h) - 2*J(z) + J(z - h))/h**2 for z in zs[1:]]
print('min J\'\' on (0, 2 ln rho]: %.4f' % min(J2))
def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
for c in [2, 3, 5, 10, 20, 47]:
    S = -D + 2*sum(math.log(ppow(n)[0])/math.sqrt(n) for n in range(2, c) if ppow(n))
    print('S(%d) = %.6f' % (c, S))
