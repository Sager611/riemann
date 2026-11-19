#!/usr/bin/env python3
"""indep_check.py -- LiRangeAudit: (I) independent frontier n* of each certified inequality (mpmath quad of the loss with exact E',
no Riemann sums), compared with the certified N (must have N <= n*); (II) no-gap sampling: the shipped certificate functions
(HP141/Li200k/verify.py = Li200k with C3 = 9.4926; DirLi/certify.py) checked at 60 log-uniform random n <= N per variant, and
at N(1+1e-6), N(1+1e-4) (expected to fail: N is at the frontier)."""
import sys, math, random, time, importlib.util
import mpmath as mp
mp.mp.dps = 30
t_start = time.time()
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
V = load('verify141', 'Riemann/Deep/HP141/Li200k/verify.py')
sys.path.insert(0, 'Riemann/Deep/DirLi'); C = load('certify', 'Riemann/Deep/DirLi/certify.py')
from fractions import Fraction
TH141 = Fraction(8471226452996444665947167626627751394950993272108499633744655722818864115887871138475393503292677558375056192430877447107716734, 22690024047233559201583164646746989164574500567295919931571742022148695950234080226695489006874468753832869742587965378849637400)
EG = mp.euler
# ---------------- (I) zeta
T0 = mp.mpf(3000175332800)
def Ez(t): return mp.mpf('0.1038') * mp.log(t) + mp.mpf('0.2573') * mp.log(mp.log(t)) + mp.mpf('9.4926') + mp.mpf(7) / 8
def dEz(t): return (mp.mpf('0.1038') + mp.mpf('0.2573') / mp.log(t)) / t
def qz(n, a, t): x = n * a / (2 * t * t); return 4 * mp.sinh(x) ** 2 + 4 * n * a * a / t ** 3
def Qz(n, af):
    f = lambda s: qz(n, af(T0 * s), T0 * s) * (mp.log(T0 * s / (2 * mp.pi)) / (2 * mp.pi) + dEz(T0 * s)) * T0
    return 2 * qz(n, af(T0), T0) * Ez(T0) + mp.quad(f, [1, 1.05, 1.2, 1.5, 2, 4, 10, 100, mp.inf])
def tz(n): return mp.cot(mp.pi / (2 * n)) / 2
def Mz(n):   # Li200k gain B-form at n1 = n2 = n (V l. 84-103)
    tl = mp.sqrt(n); U = n * 2 * mp.atan(1 / (2 * tl)); t1 = tz(n)
    main = n / mp.pi * ((mp.pi / 2 - 2 / U) * (mp.log(n / (2 * mp.pi)) - 1 / (4 * tl * tl)) - mp.pi / 2 * (1 - EG))
    return main - 8 * Ez(t1) * (U / mp.pi + 3) - (mp.mpf('0.1038') + mp.mpf('0.2573') / mp.log(t1)) * n * n / (2 * t1 * t1)
def frontier(F, lo, hi):
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    for _ in range(90):
        mid = mp.sqrt(lo * hi)
        if F(mid) > 0: lo = mid
        else: hi = mid
    return lo
const = lambda v: (lambda t: mp.mpf(v.numerator) / v.denominator)
mty = lambda t: mp.mpf(1) / 2 - 1 / (mp.mpf('5.558691') * mp.log(t))
Z = [('A a=1/2', const(Fraction(1, 2)), 708505587700793981366435840, 'A (a=1/2)'),
     ('B MTY', mty, 717734280405823147676270592, 'B (a=1/2-1/(5.558691 log t))'),
     ('C a=3/8', const(Fraction(3, 8)), 952080909647480761027133440, 'C (a=3/8, 7/8 half-plane)'),
     ('G theta141', const(TH141 - Fraction(1, 2)), 956413498941653124095410176, 'G (a=theta141-1/2, PR 141 exact, 13 Nov 2026)')]
for name, af, N, var in Z:
    ns = frontier(lambda n: Mz(n) - Qz(n, af), 10 ** 25, 10 ** 28)
    print("I zeta %-10s certified N = %.12e, independent frontier n* = %.12e, N/n* - 1 = %.3e (must be <= 0)" % (name, N, float(ns), float(mp.mpf(N) / ns - 1)), flush=True)

# ---------------- (I) DirLi, single modulus q (halved bounds, P §5B/§6)
def Ed(t, q): l = mp.log(q * (t + 2) / (2 * mp.pi)); return mp.mpf('0.22737') * l + 2 * mp.log(1 + l) - mp.mpf(1) / 2
def dEd(t, q): l = mp.log(q * (t + 2) / (2 * mp.pi)); return (mp.mpf('0.22737') + 2 / (1 + l)) / (t + 2)
def cEd(t, q): return mp.mpf('0.22737') + 2 / (1 + mp.log(q * (t + 2) / (2 * mp.pi)))
def Qd(n, a, q, t0):
    f = lambda s: qz(n, a, t0 * s) * (mp.log(q * t0 * s / (2 * mp.pi)) / mp.pi + dEd(t0 * s, q)) * t0
    return (2 * qz(n, a, t0) * Ed(t0, q) + mp.quad(f, [1, 1.05, 1.2, 1.5, 2, 4, 10, 100, mp.inf])) / 2
def Md(n, q):
    tl = mp.sqrt(n); U = n * 2 * mp.atan(1 / (2 * tl)); t1 = tz(n)
    main = n / mp.pi * ((mp.pi / 2 - 2 / U) * (mp.log(q * n / (2 * mp.pi)) - 1 / (4 * tl * tl)) - mp.pi / 2 * (1 - EG))
    return main - 4 * Ed(t1, q) * (U / mp.pi + 2) - cEd(t1, q) * n * n / (4 * t1 * t1)
D = [(3, Fraction(1, 2), 60469357461824346), (3, Fraction(3, 8), 81573243404167256), (1000, Fraction(1, 2), 413010283643), (1000, Fraction(3, 8), 559492859836)]
for q, a, N in D:
    t0 = C.t0_of(q); t0m = mp.mpf(t0.numerator) / t0.denominator; am = mp.mpf(a.numerator) / a.denominator
    ns = frontier(lambda n: Md(n, q) - Qd(n, am, q, t0m), 10 ** 4, 10 ** 20)
    print("I dirli q=%-5d a=%s certified N = %.12e, independent frontier n* = %.12e, N/n* - 1 = %.3e (must be <= 0)" % (q, a, N, float(ns), float(mp.mpf(N) / ns - 1)), flush=True)
# ---------------- (II) no-gap sampling with the shipped certificate functions
rnd = random.Random(277)
for name, af, N, var in Z:
    if var not in V.VARIANTS: continue
    bad = [n for n in sorted(int(math.exp(rnd.uniform(math.log(89), math.log(N)))) for _ in range(60)) if not (V.M_lo(n, n) - V.Q_hi(n, var) > 0)]
    beyond = [(f, bool(V.M_lo(int(N * f), int(N * f)) - V.Q_hi(int(N * f), var) > 0)) for f in (1 + 1e-6, 1 + 1e-4, 1 + 1e-2)]
    print("II verify141 %-10s 60 random n in [89, N]: failures %s; certified at N(1+e): %s" % (name, bad or 'none', beyond), flush=True)
def lowerd(n, p):
    gs = [C.gainC(n, p) if n <= 400 else None, C.gainA(n, n, p) if n >= 3 else None, C.gainB(n, n, p) if n >= 10 ** 4 else None]
    gs = [g for g in gs if g is not None]
    return max(gs, key=lambda g: float(g)) - C.loss(n, p) if gs else None
for q, a, N in D:
    p = C.P(q, q, C.t0_of(q), a, 'q=%d' % q)
    smp = sorted(set(int(math.exp(rnd.uniform(0, math.log(N)))) for _ in range(60)) | {1, 2, 3, 4, 5})
    bad = [n for n in smp if not (lowerd(n, p) is not None and lowerd(n, p) > 0)]
    beyond = [(f, bool(lowerd(int(N * f), p) > 0)) for f in (1 + 1e-6, 1 + 1e-4, 1 + 1e-2)]
    print("II certify q=%-5d a=%s %d random n in [1, N]: failures %s; certified at N(1+e): %s" % (q, a, len(smp), bad or 'none', beyond), flush=True)
print("total %.0f s" % (time.time() - t_start))

