# ceilaudit B3: detector side's one-piece limit. kappa*(xi) of Conrey's one-piece functional (zeta: xi = nu) from gl2ceil's exact
# Euler-Lagrange evaluation (Riemann/Deep/GL2Ceil/ceil.py, imported read-only; beta = Q(0)+Q(1) free, so >= SL's Q(1) = 0 class),
# then theta_c* = min over nu of max{(4 nu + 1)/3, theta_det(kappa*(nu))}.
import sys, os, numpy as np
from scipy.optimize import minimize
from mpmath import mp, mpf, cot, sqrt as msqrt, findroot
G = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'GL2Ceil'))
sys.argv = ['x', '25', '256', '0.07006']; sys.path.insert(0, G)
import ceil as cc
def kap_ode(R, lam, xi):
    cc.xi = xi; sc, mc, beta, J = cc.build(R, lam); return 1 - np.log(J)/R, beta
def argmax(xi, x0):
    r = minimize(lambda v: -kap_ode(v[0], abs(v[1]), xi)[0], x0, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-13, maxiter=400))
    return -r.fun, r.x[0], abs(r.x[1])
mp.dps = 25
cth = lambda t: 2 - t/2 - cot(t/msqrt(2))/msqrt(2)
kst = lambda t: 1 - 2/(2 - cth(t))
x0 = [8.8, 0.5]
k, R, lam = argmax(3/25, x0); print(f"check xi = 3/25: kappa* = {k:.10f} at R = {R:.4f}, lam = {lam:.4f} (frozen SL detector: 0.0844779080)", flush=True)
best = None
for nu in [0.1170, 0.1180, 0.1190, 0.1200, 0.1210]:
    k, R, lam = argmax(nu, [R, lam]); beta = kap_ode(R, lam, nu)[1]
    tm = (4*nu + 1)/3; td = float(findroot(lambda t: kst(t) - mpf(k), mpf('0.494')))
    on = max(tm, td); best = on if best is None else min(best, on)
    print(f"nu = {nu:.4f}: kappa* = {k:.10f} (R = {R:.4f}, lam = {lam:.4f}, Q(1) = beta - 1 = {beta-1:.2e}); moment theta > {tm:.6f}, detector theta > {td:.6f} -> {on:.6f}", flush=True)
print(f"one-piece detector limit on this grid: theta_c* ~ {best:.6f}")
