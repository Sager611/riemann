#!/usr/bin/env python3
"""selbaudit check 3 (DERIVATION §8): J(B) = int_1^B (B - x) e^{i rho x} d(x;u1) d(x;u2) dx, B = A*kappa, rho = Lambda/log xi in [2, 6],
d(x;u) = [Delta_2(x l1) - Delta_2(x l2)]/x exact (psi by Gauss-Legendre steps). Claim: |J| <= IBP bound (linear in B):
(1/rho)[(B-1)S1S1 + B(Dm1Dp2 + Dm2Dp1)] + (1/rho^2)[S1S1 + Dm1Dm2 + Dm1Dp2 + Dm2Dp1].
uv run --no-project --with numpy python3 ibpcheck.py"""
import math, time
import numpy as np
T0 = time.time()
gx, gw = np.polynomial.legendre.leggauss(10)
def psi_at(Y0, h, n):
    """psi(Y0 + k h), k = 0..n; psi(Y) = int_0^Y e^{-it} t^{-1/2} dt = 2 int_0^{sqrt Y} e^{-i y^2} dy."""
    out = np.zeros(n + 1, dtype=complex)
    if Y0 > 0:
        z = math.sqrt(Y0); m = 400; ed = np.linspace(0, z, m + 1)
        a = ed[:-1, None]; b = ed[1:, None]; y = (a + b)/2 + (b - a)/2*gx[None, :]
        out[0] = 2*np.sum((b - a)/2*gw[None, :]*np.exp(-1j*y**2))
    a = Y0 + h*np.arange(n)[:, None]; b = a + h; t = (a + b)/2 + (b - a)/2*gx[None, :]
    if Y0 == 0:   # first step: substitute t = y^2 on [0, h]
        z = math.sqrt(h); yy = z/2 + z/2*gx
        first = 2*np.sum(z/2*gw*np.exp(-1j*yy**2))
    inc = np.sum((b - a)/2*gw[None, :]*np.exp(-1j*t)/np.sqrt(t), axis=1)
    if Y0 == 0: inc[0] = first
    out[1:] = out[0] + np.cumsum(inc)
    return out
def D2(Y, ps):
    return -2*np.sqrt(Y)*np.exp(-1j*Y) - (1 + 2j*Y)*ps
TH = 0.15; hx = 0.005; XMAX = 1300.0
x = 1 + hx*np.arange(int((XMAX - 1)/hx) + 1); n = len(x) - 1
def dfun(u):
    l1 = 1 - u; l2 = max(TH - u, 0.0)
    p1 = psi_at(l1, l1*hx, n); Y1 = l1*x
    d = D2(Y1, p1); dp = p1.copy()
    if l2 > 0:
        p2 = psi_at(l2, l2*hx, n); Y2 = l2*x; d = d - D2(Y2, p2); dp = dp - p2
    return d/x, dp/x**2
def simpson_cum(F):
    """cumulative int_1^{x_k} F for even k (Simpson), returned at even indices."""
    s = hx/3*(F[0:-2:2] + 4*F[1:-1:2] + F[2::2]); return np.concatenate([[0], np.cumsum(s)])
print("theta=%.2f, x in [1,%g], step %g" % (TH, XMAX, hx))
Bs = [2, 5, 10, 30, 100, 300, 1000, 1290]
worst = 0.0
for (u1, u2) in [(0.0, 0.0), (0.15, 0.15), (0.0, 0.5), (0.3, 0.6), (0.05, 0.9), (0.15, 0.4)]:
    d1, dp1 = dfun(u1); d2, dp2 = dfun(u2)
    S1a, S1b = abs(d1[0]), abs(d2[0])
    Dm1 = max(np.max(abs(d1)), abs(-2*math.sqrt(math.pi)*(TH if u1 < TH else 0) + 0) ) ; Dm2 = np.max(abs(d2))
    Dm1 = np.max(abs(d1)); Dp1 = np.sum(abs(dp1[:-1]))*hx + 2*2.39/XMAX; Dp2 = np.sum(abs(dp2[:-1]))*hx + 2*2.39/XMAX
    for rho in [2.0, 3.0, 6.0]:
        F = np.exp(1j*rho*x)*d1*d2
        C0 = simpson_cum(F); C1 = simpson_cum(x*F); xe = x[::2]
        row = []
        for B in Bs:
            k = int(round((B - 1)/(2*hx))); Bk = xe[k]
            J = Bk*C0[k] - C1[k]
            bnd = ((Bk - 1)*S1a*S1b + Bk*(Dm1*Dp2 + Dm2*Dp1))/rho + (S1a*S1b + Dm1*Dm2 + Dm1*Dp2 + Dm2*Dp1)/rho**2
            crude = np.sum(np.abs((Bk - x[:2*k + 1])*d1[:2*k + 1]*d2[:2*k + 1]))*hx   # no oscillation used: ~B^2
            worst = max(worst, abs(J)/bnd)
            row.append("B=%g |J|=%.3g bnd=%.3g r=%.2f |J|/B=%.3g crude=%.3g" % (B, abs(J), bnd, abs(J)/bnd, abs(J)/Bk, crude))
        print("u=(%.2f,%.2f) rho=%.0f S1=(%.3f,%.3f) Dm=(%.3f,%.3f) Dp=(%.3f,%.3f): " % (u1, u2, rho, S1a, S1b, Dm1, Dm2, Dp1, Dp2) + " | ".join(row), flush=True)
print("max |J|/IBP-bound over all cases: %.4f  (claim: <= 1)   (%.0fs)" % (worst, time.time() - T0))
