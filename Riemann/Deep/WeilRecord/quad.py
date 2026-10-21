# pyright: reportMissingImports=false
"""Rigorous composite Gauss-Legendre quadrature on [0, R] with Bernstein-ellipse error bounds.

Panels: [0,1/2],[1/2,1],[1,2],[2,4],... (dyadic) up to length lmax, then steps of lmax to R.
For f holomorphic on the closed Bernstein ellipse E_rho of a panel of half-length h, |f| <= M there,
n-point Gauss error <= h * (64/15) * M * rho^(-2n) / (rho^2 - 1)   (Trefethen, ATAP Thm 19.3, scaled).
Ellipse of panel [c-h, c+h]: foci c +- h, semi-axes A = h(rho+1/rho)/2, B = h(rho-1/rho)/2;
min distance from the real panel to the boundary is A - h (Cauchy estimates for node rounding).
max|g| on the ellipse for holomorphic g: max over boundary arcs, each enclosed in a complex ball.
"""
import math
from flint import arb, acb, ctx

def panels(R, lmax):
    """Exact panel endpoints (arb points). R: arb (exact point) or number."""
    R = arb(R)
    ends = [arb(0), arb(1)/2, arb(1)]
    x = 1
    while 2*x <= min(float(R), lmax) + 1e-12:
        x *= 2; ends.append(arb(x))
    while float(ends[-1]) + lmax < float(R) - 1e-12:
        ends.append(ends[-1] + lmax)
    if float(ends[-1]) < float(R) - 1e-12: ends.append(R)
    else: ends[-1] = R
    return [(ends[i], ends[i+1]) for i in range(len(ends)-1)]

def ellipse_geom(lo, hi, rho):
    c = (lo + hi)/2; h = (hi - lo)/2
    rho = arb(rho)
    A = h*(rho + 1/rho)/2; B = h*(rho - 1/rho)/2
    return c, h, A, B

def ellipse_boxes(lo, hi, rho, nbox):
    """Complex balls covering the boundary of the Bernstein ellipse of [lo,hi]."""
    c, h, A, B = ellipse_geom(lo, hi, rho)
    out = []
    for k in range(nbox):
        th = 2*arb.pi()*(k + arb(1)/2)/nbox
        z = acb(c + A*th.cos(), B*th.sin())
        r = float((A*arb.pi()/nbox).upper())*(1 + 1e-12) + 1e-300   # |dz/dtheta| <= A, half-arc = pi/nbox
        out.append(acb(z.real + arb(0, r), z.imag + arb(0, r)))
    return out

def _arc_ball(c, A, B, t0, t1):
    tm = (t0 + t1)/2
    z = acb(c + A*tm.cos(), B*tm.sin())
    r = float((A*(t1 - t0)/2).upper())*(1 + 1e-12) + 1e-300   # |dz/dtheta| <= A
    return acb(z.real + arb(0, r), z.imag + arb(0, r))

def max_abs_on_ellipse(g, lo, hi, rho, nbox=64, maxdepth=14):
    """Rigorous upper bound of |g| on the boundary of the Bernstein ellipse (adaptive arc bisection)."""
    c, h, A, B = ellipse_geom(lo, hi, rho)
    m = arb(0)
    stack = [(2*arb.pi()*k/nbox, 2*arb.pi()*(k+1)/nbox, 0) for k in range(nbox)]
    while stack:
        t0, t1, d = stack.pop()
        v = abs(g(_arc_ball(c, A, B, t0, t1)))
        if v.is_finite() and float(v.rad()) <= 0.5*float(abs(v.mid())) + 1.0:
            m = max(m, arb(v.upper())); continue
        if d >= maxdepth:
            if not v.is_finite(): return None
            m = max(m, arb(v.upper())); continue
        tm = (t0 + t1)/2
        stack.append((t0, tm, d+1)); stack.append((tm, t1, d+1))
    return m

def gauss_err_factor(h, rho, n):
    rho = arb(rho)
    return h*(arb(64)/15)*rho**(-2*n)/(rho**2 - 1)

def digamma_poles_clear(lo, hi, rho):
    """True if the ellipse avoids every pole z = +-i(2k+1/2) of psi(1/4 +- i z/2)."""
    c, h, A, B = ellipse_geom(lo, hi, rho)
    if c - A > 0: return True                       # ellipse in Re z > 0
    for k in range(0, 50):
        y = arb(2*k) + arb(1)/2                       # pole at (0, +-y)
        val = (c/A)**2 + (y/B)**2                     # point (0,y) relative to centre (c,0)
        if not (val > 1): return False
    return bool(B < 100)
