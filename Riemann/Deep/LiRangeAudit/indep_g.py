#!/usr/bin/env python3
"""indep_g.py -- the G row of indep_check.py redone: a = the HP141 rational itself (it already is theta141 - 1/2 = 0.37335...);
indep_check.py subtracted 1/2 a second time (a = -0.1267), so its G line (n* = 2.9e27) is void."""
import mpmath as mp
from fractions import Fraction
mp.mp.dps = 30
A = Fraction(8471226452996444665947167626627751394950993272108499633744655722818864115887871138475393503292677558375056192430877447107716734, 22690024047233559201583164646746989164574500567295919931571742022148695950234080226695489006874468753832869742587965378849637400)
a = mp.mpf(A.numerator) / A.denominator
T0 = mp.mpf(3000175332800); EG = mp.euler
def Ez(t): return mp.mpf('0.1038') * mp.log(t) + mp.mpf('0.2573') * mp.log(mp.log(t)) + mp.mpf('9.4926') + mp.mpf(7) / 8
def dEz(t): return (mp.mpf('0.1038') + mp.mpf('0.2573') / mp.log(t)) / t
def qz(n, t): x = n * a / (2 * t * t); return 4 * mp.sinh(x) ** 2 + 4 * n * a * a / t ** 3
def Qz(n):
    f = lambda s: qz(n, T0 * s) * (mp.log(T0 * s / (2 * mp.pi)) / (2 * mp.pi) + dEz(T0 * s)) * T0
    return 2 * qz(n, T0) * Ez(T0) + mp.quad(f, [1, 1.05, 1.2, 1.5, 2, 4, 10, 100, mp.inf])
def tz(n): return mp.cot(mp.pi / (2 * n)) / 2
def Mz(n):
    tl = mp.sqrt(n); U = n * 2 * mp.atan(1 / (2 * tl)); t1 = tz(n)
    main = n / mp.pi * ((mp.pi / 2 - 2 / U) * (mp.log(n / (2 * mp.pi)) - 1 / (4 * tl * tl)) - mp.pi / 2 * (1 - EG))
    return main - 8 * Ez(t1) * (U / mp.pi + 3) - (mp.mpf('0.1038') + mp.mpf('0.2573') / mp.log(t1)) * n * n / (2 * t1 * t1)
lo, hi = mp.mpf(10) ** 26, mp.mpf(10) ** 28
for _ in range(90):
    mid = mp.sqrt(lo * hi)
    if Mz(mid) - Qz(mid) > 0: lo = mid
    else: hi = mid
N = 956413498941653124095410176
print("I zeta G a=%s certified N = %.12e, independent frontier n* = %.12e, N/n* - 1 = %.3e (must be <= 0)" % (mp.nstr(a, 12), N, float(lo), float(mp.mpf(N) / lo - 1)))
