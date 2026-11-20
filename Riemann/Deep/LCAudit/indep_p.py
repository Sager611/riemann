#!/usr/bin/env python3
"""lcaudit: independent recomputation of p(c) and 2p - 1 (c = 1e-4, 1e-3) by a different method:
K by Gauss-Legendre quadrature (not the closed form), sup's by dense grids + golden-section refinement, final values
re-evaluated with mpmath quadrature at 40 digits.  Formulas: DERIVATION §2-3 (hinge form (H)), cert.py's choice rho = 1.2036.
Run: uv run --no-project --with numpy --with mpmath python3 indep_p.py
"""
import numpy as np, mpmath as mp, math
mp.mp.dps = 40
xg, wg = np.polynomial.legendre.leggauss(400); xg = xg/2; wg = wg/2          # nodes on [-1/2, 1/2]
w0 = math.sqrt(2); Z = 2*math.sin(w0/2)/w0
f0 = np.cos(w0*xg)/Z
def Kr(u): return float(np.sum(wg*f0*np.cos(2*math.pi*u*xg)))
def Ki(y):
    y = np.atleast_1d(y)
    return (np.cosh(2*math.pi*np.outer(y, xg))*f0) @ wg
def H(y, thb):
    y = np.atleast_1d(np.asarray(y, float)); res = np.empty_like(y)
    for i in range(0, len(y), 2000):              # chunks: 2000 x 400 doubles = 6.4 MB (memory cap)
        yc = y[i:i+2000]
        res[i:i+2000] = 1 - (1 - thb)**2*4*Ki(yc)**2/(1 + Ki(2*yc))
    return res
# mpmath versions (40 digits, tanh-sinh quadrature)
mf0 = lambda x: mp.cos(mp.sqrt(2)*x)/(2*mp.sin(mp.sqrt(2)/2)/mp.sqrt(2))
mKr = lambda u: mp.quad(lambda x: mf0(x)*mp.cos(2*mp.pi*u*x), [-0.5, 0, 0.5])
mKi = lambda y: mp.quad(lambda x: mf0(x)*mp.cosh(2*mp.pi*x*y), [-0.5, 0, 0.5])
mH = lambda y, thb: 1 - (1 - thb)**2*4*mKi(y)**2/(1 + mKi(2*y))
out = []
def say(m): print(m, flush=True); out.append(m)
# C_MT: int f0^2 + iint |u-v| f0 f0 by quadrature vs closed form
I1 = mp.quad(lambda x: mf0(x)**2, [-0.5, 0.5])
I2 = mp.quad(lambda u: mp.quad(lambda v: abs(u - v)*mf0(u)*mf0(v), [-0.5, u, 0.5]), [-0.5, 0.5])
C = I1 + I2
Ccl = mp.mpf(1)/2 + mp.cot(1/mp.sqrt(2))/mp.sqrt(2)
say("C(f0) by quadrature = %s ; closed form 1/2 + cot(1/sqrt2)/sqrt2 = %s ; diff %.1e" % (mp.nstr(C, 15), mp.nstr(Ccl, 15), float(C - Ccl)))
def Bfun(rho, th=mp.mpf(4)/7):
    R = 2*mp.pi*rho
    return mp.log(1 + (1 - mp.exp(-2*R))/(mp.exp(2*th*R) - 1))/(4*mp.pi)
for c in (1e-4, 1e-3):
    tau = c
    k00, k0, k1 = mKr(tau), mKr(2*tau), mKr(4*tau)
    eps0 = 1 - k0; betaT = (6 - 16*k0/3 - 8*k1/3)/3
    lam = max(2*eps0, betaT, mp.mpf(0)); thb = mp.pi**2*tau**2/2
    assert 2*k00**2 - 1 >= 0 and 4*tau <= 1
    s = 3*float(lam)
    # rho_max: first crossing of H - s through 0 on (0.1, 3]; check H < s on (0, rho] with a dense grid
    yy = np.arange(1e-3, 3.0, 1e-5)
    hv = H(yy, float(thb)) - s
    idx = np.argmax(hv >= 0); r_lo, r_hi = yy[idx - 1], yy[idx]
    rm = mp.findroot(lambda y: mH(y, thb) - 3*lam, (r_lo + r_hi)/2)
    rho = mp.mpf('1.2036')
    assert rho < rm and np.max(hv[yy <= float(rho)]) < 0
    # alpha_rho = sup_{y > rho} (H - s)^+/(y - rho): grid on (rho, 60], then golden-section around the best grid point
    g1 = np.arange(float(rho) + 1e-5, 10.0, 1e-5); g2 = np.arange(10.0, 60.0, 1e-3)
    best = (-1.0, None)
    for g in (g1, g2):
        vals = (H(g, float(thb)) - s)/(g - float(rho))
        j = int(np.argmax(vals))
        if vals[j] > best[0]: best = (float(vals[j]), float(g[j]))
    ystar = best[1]
    gfun = lambda y: (mH(y, thb) - 3*lam)/(y - rho)
    a_, b_ = mp.mpf(ystar) - mp.mpf('2e-5'), mp.mpf(ystar) + mp.mpf('2e-5')
    gr = (mp.sqrt(5) - 1)/2
    for _ in range(60):
        c1_, c2_ = b_ - gr*(b_ - a_), a_ + gr*(b_ - a_)
        if gfun(c1_) > gfun(c2_): b_ = c2_
        else: a_ = c1_
    ys = (a_ + b_)/2; alpha = gfun(ys)
    tail = (1 - 3*lam)/(60 - rho)
    assert alpha > tail
    Bv = Bfun(rho)
    p = (2 - lam - Ccl - alpha*Bv)/(1 - lam)
    say("c = %g: lambda = %s (betaT = %s), thb = %.3e, rho_max = %s, rho = 1.2036, y* = %s, alpha_rho = %s, B(rho) = %s" %
        (c, mp.nstr(lam, 6), mp.nstr(betaT, 4), float(thb), mp.nstr(rm, 8), mp.nstr(ys, 8), mp.nstr(alpha, 10), mp.nstr(Bv, 8)))
    say("        p = %s, 2p - 1 = %s" % (mp.nstr(p, 12), mp.nstr(2*p - 1, 12)))
    # rho scan: is rho = 1.2036 near-optimal?
    best_r = None
    for rr in np.arange(1.10, float(rm) - 1e-4, 0.005):
        g = np.arange(rr + 1e-4, 10.0, 2e-4)
        a_r = float(np.max((H(g, float(thb)) - s)/(g - rr)))
        pr = (2 - float(lam) - float(Ccl) - a_r*float(Bfun(mp.mpf(rr))))/(1 - float(lam))
        if best_r is None or pr > best_r[0]: best_r = (pr, rr, a_r)
    say("        rho scan (grid alpha): best p = %.9f at rho = %.3f (alpha %.6f)" % best_r)
open('logs/indep_p.log', 'w').write('\n'.join(out) + '\n')
