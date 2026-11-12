"""Twist check: e^{ic(r2-r1)/2} * U_L conj(U_R) = 2cos(c(X1-X2-X3+X4)/2) - 2cos(c(X1-X2+X3-X4)/2 + c sigma) is REAL, so with
f = e^{icU/2} g the complex problem (orientation L) is the REAL problem with this kernel.  Check: Im of the twisted Taylor arrays ~ 0
and the real twisted pencil reproduces the complex pencil value. usage: twcheck.py r kappa N"""
import sys, time
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
sys.path.insert(0, HERE)
import gcore as G
import gcore_r as GR
import ctest as CT
r = int(sys.argv[1]); kap = float(sys.argv[2]); N = int(sys.argv[3]); c = kap * np.pi
orig = G.rho_parts
def rho_tw(X1, X2, X3, X4, c):
    out = orig(X1, X2, X3, X4, c)
    tw = np.exp(1j * c * (X3 + X4 - X1 - X2) / 2)
    return [o * tw for o in out]
a, b, cc = CT.arrays(c, r, 12)
G.rho_parts = rho_tw
at, bt, cct = CT.arrays(c, r, 12)
for nm, h in (('a', at), ('b', bt), ('cc', cct)):
    hs = h[:6, :6]
    print('twisted %s: max|Im| (i,j<6) = %.2e, max|h| = %.2e, max|h-h^T| = %.2e' % (nm, np.abs(hs.imag).max(), np.abs(hs).max(), np.abs(hs - hs.T).max()))
print('raw a[0,1] = %s ; twisted a[0,1] = %s' % (a[0, 1], at[0, 1]))
t0 = time.time()
M0, MS, MSS, _ = CT.mats(c, 0.25, r, N, 36, False)   # real part of twisted arrays (G.rho_parts patched)
lt, vt = CT.best(M0, MS, MSS)
G.rho_parts = orig
M0, MS, MSS, _ = CT.mats(c, 0.25, r, N, 36, True)
lc, vc = CT.best(M0, MS, MSS)
print('r=%d kappa=%.4f N=%d: twisted-real pencil %.7f (v=%.4f) | complex pencil %.7f (v=%.4f) [%.0fs]' % (r, kap, N, lt / kap ** 2, vt, lc / kap ** 2, vc, time.time() - t0))
