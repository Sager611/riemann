#!/usr/bin/env python3
"""identities.py -- symbolic / exact / high-precision checks of the formulas used in PROOF.md (Stage A).
Prints 'ok'/'FAIL' lines and 'IDENTITIES ALL OK' at the end."""
import sys
import sympy as sp
import mpmath as mp
from fractions import Fraction
from math import comb, factorial
FAIL = []
def check(c, msg):
    print(('ok   ' if c else 'FAIL ') + msg, flush=True)
    if not c: FAIL.append(msg)
z, p, r, ph = sp.symbols('z p r phi', positive=True)
# (1) Liouville transformation of w'' = (nu^2 f + g) w,  f = (1+z^2)/z^2, g = -1/(4 z^2)
f = (1 + z**2) / z**2; g = -1 / (4 * z**2)
fm14 = sp.sqrt(z) * (1 + z**2)**sp.Rational(-1, 4); f34 = z**sp.Rational(-3, 2) * (1 + z**2)**sp.Rational(3, 4)
psi = sp.simplify(g / f - sp.diff(fm14, z, 2) / f34)
check(sp.simplify(psi - z**2 * (4 - z**2) / (4 * (1 + z**2)**3)) == 0, 'psi = z^2(4-z^2)/(4(1+z^2)^3)')
U1 = (3 * p - 5 * p**3) / 24
dU1 = sp.diff(U1.subs(p, (1 + z**2)**sp.Rational(-1, 2)), z)
check(sp.simplify(psi * sp.sqrt(1 + z**2) / z - 2 * dU1) == 0, 'psi d xi = 2 dU1(p),  d xi/dz = sqrt(1+z^2)/z')
eta = sp.sqrt(1 + z**2) + sp.log(z / (1 + sp.sqrt(1 + z**2)))
check(sp.simplify(sp.diff(eta, z) - sp.sqrt(1 + z**2) / z) == 0, "eta' = sqrt(1+z^2)/z")
# (2) D = z d/dz identities for Psi_e = sqrt(1+z^2) - log(1+sqrt(1+z^2)) and P = log(1+z^2)/4
Pe = sp.sqrt(1 + z**2) - sp.log(1 + sp.sqrt(1 + z**2)); s = sp.sqrt(1 + z**2)
D = lambda e: sp.simplify(z * sp.diff(e, z))
d1 = D(Pe); d2 = D(d1); d3 = D(d2); d4 = D(d3)
check(sp.simplify(d1 - z**2 / (s + 1)) == 0 and sp.simplify(d1 - (s - 1)) == 0, 'D Psi_e = s - 1 = z^2/(s+1)')
check(sp.simplify(d2 - z**2 / s) == 0, 'D^2 Psi_e = z^2/s')
check(sp.simplify(d3 - z**2 * (2 + z**2) / s**3) == 0, 'D^3 Psi_e = z^2(2+z^2)/s^3')
check(sp.simplify(d4 - z**2 * (1 / s + 3 / s**5)) == 0, 'D^4 Psi_e = z^2(1/s + 3/s^5)')
P = sp.log(1 + z**2) / 4; u = 1 / (1 + z**2); v = z**2 * u
q1 = D(P); q2 = D(q1); q3 = D(q2); q4 = D(q3)
check(sp.simplify(q1 - v / 2) == 0 and sp.simplify(q2 - u * v) == 0, 'D P = v/2, D^2 P = u v  (u = 1/(1+z^2), v = z^2 u)')
check(sp.simplify(q3 - (-2 * u * v * (1 - 2 * u))) == 0, 'D^3 P = -2uv(1-2u)')
check(sp.simplify(q4 - 4 * u * (1 - 7 * u + 12 * u**2 - 6 * u**3)) == 0 and sp.simplify(q4 - 4 * u * v * (1 - 6 * v + 6 * v**2)) == 0, 'D^4 P = 4u(1-7u+12u^2-6u^3) = 4uv(1-6v+6v^2)')
# (3) d/dphi Re Psi_e(r e^{i phi}) = -Im sqrt(1 + r^2 e^{2 i phi})   (numerical check at sample points)
mp.mp.dps = 40
def PsieN(w): s = mp.sqrt(1 + w * w); return s - mp.log(1 + s)
dev = 0
for (rr, pp) in [(0.3, 0.2), (0.7, 1.0), (1.2, 0.5), (1.55, 0.4), (0.8, 1.5)]:
    rr, pp = mp.mpf(rr), mp.mpf(pp)
    der = mp.diff(lambda t: mp.re(PsieN(rr * mp.expj(t))), pp)
    dev = max(dev, abs(der + mp.im(mp.sqrt(1 + rr**2 * mp.expj(2 * pp)))))
check(dev < 1e-15, 'd/dphi Re Psi_e(r e^{i phi}) = -Im sqrt(1+r^2 e^{2 i phi})  (max dev %.1e)' % float(dev))
# (4) Pade remainder identity e^z q(z) = p(z) - z^{2m+1} E(z) (odd m), exact power series to degree 4m
def ser_check(m):
    L = 4 * m + 6
    qc = [Fraction((-1)**j * comb(m, j) * factorial(2 * m - j), factorial(2 * m)) for j in range(m + 1)] + [Fraction(0)] * L
    pc = [abs(x) for x in qc]
    ez = [Fraction(1, factorial(k)) for k in range(L)]
    lhs = [sum(ez[k - j] * qc[j] for j in range(0, min(k, m) + 1)) for k in range(L)]
    Ec = [Fraction(factorial(m) * factorial(m + j), factorial(j) * factorial(2 * m) * factorial(2 * m + j + 1)) for j in range(L)]
    rhs = [pc[k] - (Ec[k - 2 * m - 1] if k >= 2 * m + 1 else 0) for k in range(L)]
    return lhs == rhs
check(all(ser_check(m) for m in (1, 3, 5, 7, 9, 11)), 'e^z q = p - z^{2m+1} E  (exact series, m = 1..11 odd)')
# (5) Bessel forms (I1), (I2) at sample points: p(2 nu zeta) and X = z^{2m+1}E/p vs I_nu, K_nu
mp.mp.dps = 60
def I12(m):
    nu = mp.mpf(m) + mp.mpf(1) / 2
    pc = [mp.mpf(comb(m, j) * factorial(2 * m - j)) / factorial(2 * m) for j in range(m + 1)]
    Ec = [mp.mpf(factorial(m) * factorial(m + j)) / (factorial(j) * factorial(2 * m) * factorial(2 * m + j + 1)) for j in range(400)]
    poly = lambda c, w: mp.polyval(c[::-1], w)
    th0 = mp.mpf(factorial(2 * m)) / (factorial(m) * 2**m)
    worst = 0
    for zeta in [mp.mpf('0.3'), mp.mpc('0.5', '0.4'), mp.mpc('0.2', '0.8'), mp.mpf('1.3')]:
        x = nu * zeta; zz = 2 * x
        pv = poly(pc, zz); Ev = poly(Ec, zz)
        pB = mp.sqrt(2 / mp.pi) * x**nu * mp.exp(x) * mp.besselk(nu, x) / th0
        XB = mp.pi * mp.besseli(nu, x) / mp.besselk(nu, x)
        worst = max(worst, abs(pv / pB - 1), abs(zz**(2 * m + 1) * Ev / pv / XB - 1))
    return worst
w = max(I12(m) for m in (3, 7, 21))
check(w < mp.mpf(10)**-30, 'p(z) = sqrt(2/pi) x^nu e^x K_nu(x)/theta_m(0) and z^{2m+1}E/p = pi I_nu(x)/K_nu(x), x=z/2 (m=3,7,21; rel.err %.1e)' % float(w))
# (6) kappa = sqrt(2 pi nu) nu^nu e^-nu / Gamma(nu+1) in (exp(-1/(12 nu)), 1)  (Stirling; check)
okk = True
for nu in [mp.mpf(203) + mp.mpf(1)/2, mp.mpf(1001) + mp.mpf(1)/2, mp.mpf(10)**6 + mp.mpf(1)/2]:
    kap = mp.sqrt(2 * mp.pi * nu) * nu**nu * mp.exp(-nu) / mp.gamma(nu + 1)
    okk = okk and (mp.exp(-1 / (12 * nu)) < kap < 1)
check(okk, 'exp(-1/(12 nu)) < kappa < 1 (Stirling; spot check)')
# (7) saddle formula: r*(a,t) solves r + a(sqrt(1+r^2)-1) = t
aa, tt = sp.symbols('a t', positive=True)
rs = (-(tt + aa) + aa * sp.sqrt((tt + aa)**2 - aa**2 + 1)) / (aa**2 - 1)
chk = [float((rs + aa * (sp.sqrt(1 + rs**2) - 1) - tt).subs({aa: sp.Integer(A_), tt: sp.Rational(T_)}).evalf(40)) for (A_, T_) in [(3, '1/10'), (5, 2), (41, 9), (1001, 300)]]
check(max(abs(c) for c in chk) < 1e-20, 'saddle formula r*(a,t) (a >= 3); a = 1: r* = t(t+2)/(2(t+1))')
# (8) Lemma 1 numerically: |h_K|, |h_I| (from mpmath Bessel functions) <= the lg_box bounds of fh2lib
import fh2lib
from flint import arb
mp.mp.dps = 50
worst = 0
for nu_ in (mp.mpf('5.5'), mp.mpf('21.5'), mp.mpf('203.5')):
    kap = mp.sqrt(2 * mp.pi * nu_) * nu_**nu_ * mp.exp(-nu_) / mp.gamma(nu_ + 1)
    for (rr, th) in [(0.3, 0.0), (0.64, 0.675), (0.82, 1.33), (1.3, 0.0), (1.17, 0.54), (0.84, 1.51), (1.6, 0.3)]:
        zeta = mp.mpf(rr) * mp.expj(mp.mpf(th)); s_ = mp.sqrt(1 + zeta**2)
        eta_ = s_ - mp.log(1 + s_) + mp.log(zeta)
        hK = mp.besselk(nu_, nu_ * zeta) / (mp.sqrt(mp.pi / (2 * nu_)) * mp.exp(-nu_ * eta_) * (1 + zeta**2)**mp.mpf(-0.25)) - 1
        hI = mp.besseli(nu_, nu_ * zeta) / (kap * (2 * mp.pi * nu_)**mp.mpf(-0.5) * mp.exp(nu_ * eta_) * (1 + zeta**2)**mp.mpf(-0.25)) - 1
        lg = fh2lib.lg_box(arb(rr), arb(rr), arb(th), arb(th), arb(str(nu_)))
        if lg is None: continue
        eI, eK = float(lg[0].upper()), float(lg[1].upper())
        worst = max(worst, float(abs(hK)) / eK, float(abs(hI)) / eI)
check(worst <= 1, 'Lemma 1 (LG error bounds) holds at 21 sample points, nu = 5.5, 21.5, 203.5 (max |h|/bound = %.3f)' % worst)
print('IDENTITIES ALL OK' if not FAIL else 'IDENTITIES FAILURES: %d' % len(FAIL))
sys.exit(1 if FAIL else 0)
