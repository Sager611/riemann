import time, flint, mpmath
from flint import acb, arb, ctx
ctx.prec = 64
print([m for m in dir(acb) if 'zeta' in m or 'dirichlet' in m or 'hurwitz' in m][:40])
print([m for m in dir(flint) if 'series' in m or 'dirichlet' in m])
s = acb(0.37, 10000.5)
t0=time.time()
for k in range(200): z = acb(0.37, 10000.5 + 0.01*k).zeta()
print('acb.zeta per call', (time.time()-t0)/200, z)
try:
    t0=time.time()
    for k in range(50):
        z = flint.acb_series([acb(0.37, 10000.5+0.01*k), 1], prec=2).zeta()
    print('acb_series zeta per call', (time.time()-t0)/50, z)
except Exception as e: print('acb_series zeta err', e)
mpmath.mp.dps = 15
t0=time.time()
for k in range(50): w = mpmath.zeta(mpmath.mpc(0.37, 10000.5+0.01*k), 1, 1)
print('mpmath zeta deriv per call', (time.time()-t0)/50, w)
t0=time.time()
for k in range(50): w = mpmath.zeta(mpmath.mpc(0.37, 10000.5+0.01*k))
print('mpmath zeta per call', (time.time()-t0)/50, w)
try:
    chi = flint.dirichlet_char(5, 2)
    t0=time.time()
    for k in range(20): v = acb(0.37, 10000.5+0.01*k).dirichlet_l(chi)
    print('acb dirichlet_l per call', (time.time()-t0)/20, v)
except Exception as e: print('dirichlet err', e)
