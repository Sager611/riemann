#!/usr/bin/env python3
"""lincomb3 numerics step 3: exact zero counts of F = sum_j b_j e^{i theta_j} L(s, chi_j) (q = 7, odd chi_3, chi_5, chi_6 in Conrey labels)
by the argument principle on the right half of the rectangle [1-s1, s1] x [t0, T] (Lambda_F(1-conj s) = conj Lambda_F(s)):
N_F(t0, T) = (1/pi) [vartheta(T) - vartheta(t0) + Delta arg F along 1/2+it0 -> s1+it0 -> s1+iT -> 1/2+iT].
The phases e^{i theta_j} are read off Arb's Hardy Z (Z_j = e^{i(vartheta+theta_j)} L_j), so G = e^{i vartheta} F = sum b_j Z_j as in analyse.py.
Run: uv run --no-project --with python-flint --with mpmath --with numpy python3 argcount.py
"""
import math, numpy as np, mpmath as mp
from flint import dirichlet_char, acb, ctx
ctx.prec = 60
q, a = 7, 1
ns = (3, 5, 6)
chis = [dirichlet_char(q, n) for n in ns]
def vartheta(t):  # arg of (q/pi)^{(1/2+a+it)/2} Gamma((1/2+a+it)/2), continuous in t
    return float(t / 2 * math.log(q / math.pi) + mp.im(mp.loggamma((0.5 + a + 1j * t) / 2)))
def Lval(c, s):
    v = c.l_function(acb(s.real, s.imag)); return complex(float(v.real.mid()), float(v.imag.mid()))
# phases
ph = []
for c in chis:
    vals = []
    for tt in (100.3, 777.7):
        Zt = float(c.hardy_z(tt).real.mid()); Lt = Lval(c, complex(0.5, tt))
        e = Zt / (complex(math.cos(vartheta(tt)), math.sin(vartheta(tt))) * Lt)
        vals.append(e)
    assert abs(abs(vals[0]) - 1) < 1e-6 and abs(vals[0] - vals[1]) < 1e-6, vals
    ph.append(vals[0])
print('phases e^{i theta_j}:', ['%.6f%+.6fi' % (p.real, p.imag) for p in ph])
NMAX = 20000
nn = np.arange(1, NMAX + 1)
chiv = np.array([[complex(float(c.chi_exponent(int(m)) is not None) and 0) for m in []] for c in chis]) if False else None
def chi_values(c):
    out = np.zeros(NMAX, complex)
    for m in range(1, q + 1):
        if math.gcd(m, q) != 1: continue
        e = c.chi_exponent(m) if hasattr(c, 'chi_exponent') else None
    return out
# character values via L-function-free route: chi(m) from Arb is not exposed directly; use exponents of the Conrey label
# python-flint dirichlet_char.chi_exponent(m) gives k with chi(m) = e(k / order)
from flint import dirichlet_group
EXPO = int(dirichlet_group(q).exponent())     # chi(m) = e(chi_exponent(m) / EXPO)
cv = []
for c in chis:
    v = np.zeros(q, complex)
    for m in range(1, q):
        k = int(c.chi_exponent(m)); v[m] = complex(math.cos(2 * math.pi * k / EXPO), math.sin(2 * math.pi * k / EXPO))
    cv.append(v)
def Fdir(bv, s):  # Dirichlet series, valid for Re s >= 3
    coef = sum(bv[j] * ph[j] * cv[j][nn % q] for j in range(3))
    return complex(np.sum(coef * np.exp(-s * np.log(nn))))
def Fval(bv, s):
    if s.real >= 3.0: return Fdir(bv, s)
    return sum(bv[j] * ph[j] * Lval(chis[j], s) for j in range(3))
# check the series against Arb at s = 3 + 10i
s = complex(3.0, 10.0)
d = abs(Fdir((1, 1, 1), s) - sum(ph[j] * Lval(chis[j], s) for j in range(3)))
print('Dirichlet series check at 3+10i: |diff| = %.2e' % d)
def darg_path(bv, pts):
    v = np.array([Fval(bv, p) for p in pts])
    ang = np.unwrap(np.angle(v))
    step = np.max(np.abs(np.diff(ang)))
    return ang[-1] - ang[0], step, np.min(np.abs(v))
t0 = 5.0
NF = {}
for bv in ((1.0, 1.0, 1.0), (1.0, -2.0, 0.5), (3.0, 1.0, -1.0)):
    a1 = sum(bv[j] * ph[j] for j in range(3))
    # sigma1: first coefficient dominates the tail
    s1 = 3.0
    while abs(a1) <= 1.5 * sum(abs(x) for x in bv) * (float(mp.zeta(s1)) - 1): s1 += 1.0
    sig = np.r_[np.arange(0.5, 1.6, 0.004), np.arange(1.6, s1 + 1e-9, 0.02)]
    bot, st0, m0 = darg_path(bv, [complex(x, t0) for x in sig])
    print('\n== b = %s: |a_F(1)| = %.4f, sigma1 = %.0f; bottom edge: darg %.4f (max step %.3f)' % (bv, abs(a1), s1, bot, st0))
    Fr0 = Fval(bv, complex(s1, t0))
    for T in (600.0, 1200.0, 2400.0):
        T = T + 0.0137                     # avoid grid coincidences
        FrT = Fval(bv, complex(s1, T))
        right = math.atan2((FrT / a1).imag, (FrT / a1).real) - math.atan2((Fr0 / a1).imag, (Fr0 / a1).real)
        top, st, mn = darg_path(bv, [complex(x, T) for x in sig[::-1]])
        tot = vartheta(T) - vartheta(t0) + bot + right + top
        print('T = %7.2f: N_F(t0,T) = %.4f  (top edge max step %.3f, min |F| %.2e)' % (T, tot / math.pi, st, mn))
        NF[(bv, round(T))] = round(tot / math.pi)
    # proportion on the line: sign changes of G = sum b_j Z_j on the zgrid grid (Z_q7_n*.npy) over the exact count
    try:
        Zs = np.array([np.load('Z_q7_n%d.npy' % n) for n in ns]); tg = 5.0 + 0.03 * np.arange(Zs.shape[1])
        G = (np.array(bv)[:, None] * Zs).sum(0); sg = np.sign(G); ix = np.nonzero(sg[:-1] * sg[1:] < 0)[0]
        zc = tg[ix]
        n1, n2 = int(np.sum(zc <= 1200.0137)), int(np.sum(zc <= 2400.0137))
        a1, a2 = NF[(bv, 1200)], NF[(bv, 2400)]
        print('on-line zeros (odd multiplicity) / all zeros: (5, 2400]: %d/%d = %.4f; (1200, 2400]: %d/%d = %.4f' % (n2, a2, n2 / a2, n2 - n1, a2 - a1, (n2 - n1) / (a2 - a1)))
    except FileNotFoundError:
        print('(Z_q7_n*.npy not found: run zgrid.py first for the proportions)')
