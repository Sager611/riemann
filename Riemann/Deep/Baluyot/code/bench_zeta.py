# benchmark zeta evaluation speed at large height (flint acb vs mpmath)
import time
import flint
from flint import acb, arb, ctx
import mpmath

ctx.prec = 64
for t0 in [1e3, 1e4, 1e5, 2e5]:
    s = acb(0.6, t0)
    n = 200
    a = time.time()
    for k in range(n):
        z = acb(0.6, t0 + 0.1 * k).zeta()
    b = time.time()
    print("flint t=%g: %.1f us/eval, value %s" % (t0, (b - a) / n * 1e6, z))
    mpmath.mp.dps = 15
    a = time.time()
    m = 50
    for k in range(m):
        w = mpmath.zeta(mpmath.mpc(0.6, t0 + 0.1 * k))
    b = time.time()
    print("mpmath t=%g: %.1f us/eval, value %s" % (t0, (b - a) / m * 1e6, w))
