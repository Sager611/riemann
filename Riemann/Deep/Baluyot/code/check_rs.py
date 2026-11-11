# Accuracy of zeta_rs against Arb (python-flint) at random points
import numpy as np
from flint import acb, ctx
from zeta_rs import zeta_rs

ctx.prec = 64
rng = np.random.default_rng(1)
for T in [1e3, 1e4, 1e5]:
    for sigma in [0.5, 0.55, 0.6, 0.7, 0.8]:
        t = T + T * rng.random(40)
        z = zeta_rs(sigma, t)
        ref = np.array([complex(acb(sigma, float(x)).zeta()) for x in t])
        err = np.abs(z - ref)
        print("T=%g sigma=%.2f max|err|=%.2e median|err|=%.2e  (t^-(sigma/2+1/2) at T: %.1e)"
              % (T, sigma, err.max(), np.median(err), T ** (-(sigma / 2 + 0.5))))
