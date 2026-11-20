# rademacher.py -- numerical sanity check of the input (I2) as quoted (KLN TeX l.396-403): ratio |zeta(s)| / RHS on t-grids.
import math, time
from flint import acb, arb, ctx
import mpmath
ctx.prec = 64
t0 = time.time()
for eta in [0.1, 0.05, 0.5]:
    z1e = float(mpmath.zeta(1 + eta))
    for sig in [-eta, 0.5, 0.85, 0.92, 1.0 + eta]:
        mx = (0, None)
        for k in range(0, 20001):
            t = 0.013 + 0.5*k
            s = complex(sig, t); z = abs(complex(acb(sig, t).zeta()))
            rhs = 3*abs(1 + s)/abs(1 - s)*(abs(1 + s)/(2*math.pi))**((1 + eta - sig)/2)*z1e
            if z/rhs > mx[0]: mx = (z/rhs, t)
        print('eta', eta, 'sigma', sig, 'max |zeta|/RHS on t in [0.013, 1e4] (step 0.5) =', '%.4f' % mx[0], 'at t =', round(mx[1], 3), round(time.time() - t0), 's', flush=True)
