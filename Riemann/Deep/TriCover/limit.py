#!/usr/bin/env python3
"""Rigorous enclosure for the continuous chain model of arXiv 2609.06558 (Section 7) at the
paper's displayed (c, d).  Needs mpmath (interval arithmetic, mpmath.iv); about 11 s.
Usage: python3 limit.py

Certifies, with outward-rounded interval arithmetic at 60 digits:
 (1) a box B = tau0 +- r, s0 +- r containing exactly one zero p = (tau, s) of
       F1 = log(Rf/R0) - D s,  F2 = log(Yf/Y0) + kappa s   (Section 7 endpoint equations)
     (Krawczyk test K(B) inside int(B); this also gives det DF(p) != 0);
 (2) the strict domain inequalities of Section 5 (used by Theorem 5.3) and kappa < v on B;
 (3) g(t) = kappa R0 e^{Dt} + D Y0 e^{-kappa t} - 1 < 0 on [0, s], so E = D kappa / g > 0 is smooth;
 (4) an upper bound for S_cont = 3 (1 - c - d + int_0^s E dt) (trapezoid rule plus the
     error term (s h^2 / 12) max|E''|, both enclosed).
NOTE.md, section 4, proves that this gives lambda_simplex >= 1/S_cont(p) > 1.0633083237.
Printed bounds are rounded outward (upper bounds up, lower bounds down), so they stay valid as printed.
"""
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from mpmath import mp, iv, mpf

mp.dps = 60
iv.dps = 60

C_DEC, D_DEC = 0.3104341279739202, 0.41586782506519804   # paper, Section 7; used as these exact binary64 rationals
c = iv.mpf(C_DEC); d = iv.mpf(D_DEC)
assert c.a == c.b and d.a == d.b

class Dual:
    """value + derivative (w.r.t. one variable), both intervals."""
    def __init__(s, v, e=0):
        s.v = v if isinstance(v, type(iv.mpf(0))) else iv.mpf(v)
        s.e = e if isinstance(e, type(iv.mpf(0))) else iv.mpf(e)
    def _w(o): return o if isinstance(o, Dual) else Dual(o)
    def __add__(s, o): o = Dual._w(o); return Dual(s.v + o.v, s.e + o.e)
    __radd__ = __add__
    def __sub__(s, o): o = Dual._w(o); return Dual(s.v - o.v, s.e - o.e)
    def __rsub__(s, o): return Dual._w(o) - s
    def __mul__(s, o): o = Dual._w(o); return Dual(s.v * o.v, s.e * o.v + s.v * o.e)
    __rmul__ = __mul__
    def __truediv__(s, o): o = Dual._w(o); return Dual(s.v / o.v, (s.e * o.v - s.v * o.e) / (o.v * o.v))
    def __rtruediv__(s, o): return Dual._w(o) / s
    def __neg__(s): return Dual(-s.v, -s.e)
def dlog(x): return Dual(iv.log(x.v), x.e / x.v)
def dexp(x): ev = iv.exp(x.v); return Dual(ev, x.e * ev)

def params(tau):
    c_, d_ = (Dual(c), Dual(d)) if isinstance(tau, Dual) else (c, d)
    return _params(tau, c_, d_)

def _params(tau, c, d):
    den = 1 + tau + tau * tau
    v = (1 - d + d * tau) / den; u = tau * v; kap = 1 - u - v; D = v - d
    E0 = d * (kap + D) - kap * c; R0 = (v - c) / E0; Y0 = d / E0
    Ef = kap * (v - u) + D * (v - kap); Rf = (v - u) / Ef; Yf = v / Ef
    return u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf

def F(tau, s):
    u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf = params(tau)
    return dlog(Rf / R0) - D * s, dlog(Yf / Y0) + kap * s

def F_and_J(tau, s):
    # tau, s intervals; returns F (intervals) and Jacobian (intervals)
    f1t, f2t = F(Dual(tau, 1), Dual(s, 0))
    f1s, f2s = F(Dual(tau, 0), Dual(s, 1))
    return (f1t.v, f2t.v), ((f1t.e, f1s.e), (f2t.e, f2s.e))

def lo(x): return mpf(x.a)
def hi(x): return mpf(x.b)
def mid(x): return (lo(x) + hi(x)) / 2
def dec(x, n, up):
    """the binary number x as a decimal with n significant digits, rounded up (up=True) or down"""
    man, e = mpf(x).man_exp
    q = Fraction(man) * Fraction(2) ** e
    with localcontext() as ctx:
        ctx.prec, ctx.rounding = n, (ROUND_CEILING if up else ROUND_FLOOR)
        return str(Decimal(q.numerator) / Decimal(q.denominator))

# float Newton (mp) for the centre
tau0, s0 = mpf('0.6737630444333738'), mpf('1.1213014301809736')
for _ in range(60):
    (f1, f2), J = F_and_J(iv.mpf(tau0), iv.mpf(s0))
    f1, f2 = mid(f1), mid(f2)
    a, b, cc, dd = [mid(J[i][j]) for i in range(2) for j in range(2)]
    det = a * dd - b * cc
    dt = (dd * f1 - b * f2) / det; ds = (-cc * f1 + a * f2) / det
    tau0 -= dt; s0 -= ds
    if abs(dt) + abs(ds) < mpf(10) ** -55: break
# Krawczyk
r = mpf(10) ** -30
Bt = iv.mpf([tau0 - r, tau0 + r]); Bs = iv.mpf([s0 - r, s0 + r])
(f1, f2), _ = F_and_J(iv.mpf(tau0), iv.mpf(s0))
_, JB = F_and_J(Bt, Bs)
a, b, cc, dd = [mid(JB[i][j]) for i in range(2) for j in range(2)]
det = a * dd - b * cc
Cm = ((dd / det, -b / det), (-cc / det, a / det))
Cm = tuple(tuple(iv.mpf(x) for x in row) for row in Cm)
M = [[(1 if i == j else 0) - (Cm[i][0] * JB[0][j] + Cm[i][1] * JB[1][j]) for j in range(2)] for i in range(2)]
dx = (Bt - iv.mpf(tau0), Bs - iv.mpf(s0))
K = [ iv.mpf((tau0, s0)[i]) - (Cm[i][0] * f1 + Cm[i][1] * f2) + M[i][0] * dx[0] + M[i][1] * dx[1] for i in range(2)]
inside = (lo(K[0]) > lo(Bt) and hi(K[0]) < hi(Bt) and lo(K[1]) > lo(Bs) and hi(K[1]) < hi(Bs))
normM = max(hi(abs(M[i][0])) + hi(abs(M[i][1])) for i in range(2))
print('Krawczyk K(B) in int(B):', inside, ' ||I - C J(B)|| <=', dec(normM, 3, up=True))
assert inside and normM < 1

# domain inequalities on B
u, v, kap, D, E0, R0, Y0, Ef, Rf, Yf = params(Bt)
checks = {'tau in (0,1)': (Bt, 1 - Bt), 'u > 0': (u,), 'v - u > 0': (v - u,), 'd - v > 0': (d - v,),
          '1 - u - d > 0': (1 - u - d,), 'u + v - d > 0': (u + v - d,), 'c > 0': (c,), 'v - c > 0': (v - c,),
          'A0 = 1 - u - c - d > 0': (1 - u - c - d,), 'E0 > 0': (E0,), 'v - kappa > 0': (v - kap,), 's > 0': (Bs,)}
for name, vals in checks.items():
    ok = all(lo(x) > 0 for x in vals)
    print('  %-24s %s  (min %s)' % (name, ok, dec(min(lo(x) for x in vals), 4, up=False)))
    assert ok

# integral bound
sh = hi(Bs)
def g_parts(t):
    eD = iv.exp(D * t); ek = iv.exp(-kap * t)
    g = kap * R0 * eD + D * Y0 * ek - 1
    g1 = kap * D * R0 * eD - kap * D * Y0 * ek
    g2 = kap * D * D * R0 * eD + kap * kap * D * Y0 * ek
    return g, g1, g2
NP = 200
E2max = mpf(0)
for j in range(NP):
    t = iv.mpf([j, j + 1]) * sh / NP   # outward enclosure of [s j/NP, s (j+1)/NP]: the pieces cover all of [0, s]
    g, g1, g2 = g_parts(t)
    assert hi(g) < 0, 'denominator not negative'
    E2 = -D * kap * (g2 * g - 2 * g1 * g1) / (g * g * g)
    E2max = max(E2max, hi(abs(E2)))
N = 20000
h = iv.mpf(sh) / N
T = iv.mpf(0)
for j in range(N + 1):
    t = h * j
    g = g_parts(t)[0]
    w = iv.mpf(0.5) if j in (0, N) else iv.mpf(1)
    T += w * D * kap / g
T = T * h
# trapezoid sum plus the error term (s h^2/12) max|E''|, added in outward-rounded interval arithmetic
Iup = hi(T + iv.mpf(sh) * h * h / 12 * iv.mpf(E2max))
SupI = 3 * (1 - c - d + iv.mpf(Iup))
Sup = hi(SupI)
lam_low = iv.mpf(1) / SupI
print('int_0^s E dt <=', dec(Iup, 15, up=True), '  max|E2| <=', dec(E2max, 4, up=True))
print('S_cont <=', dec(Sup, 16, up=True))
print('lambda_simplex >= 1/S_cont >=', dec(lo(lam_low), 16, up=False))
TARGET = '1.0633083237'
assert lo(lam_low) > hi(iv.mpf(TARGET))   # iv.mpf(TARGET) encloses the decimal; compare with its upper end
print('certified: lambda_simplex >= ' + TARGET)

