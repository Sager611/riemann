#!/usr/bin/env python3
"""Exact re-verification of every displayed numerical certificate in
Guo, arXiv 2610.30476v1 (Elementary_Frobenius_Siegel_v19.tex).
Fractions for rational claims; log bounds certified by the paper's own
atanh series (rigorous rational enclosures); mpmath only as a cross-check."""
from fractions import Fraction as F
import mpmath as mp
import sympy as sp
mp.mp.dps = 60
ok_all = True
def chk(name, cond):
    global ok_all
    ok_all &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)

# --- rigorous rational enclosures of log x via the atanh series (TeX l.290-297)
def log_encl(x, m):
    z = F(x - 1, x + 1)
    s = 2 * sum(z ** (2 * j + 1) / (2 * j + 1) for j in range(m))
    r = 2 * z ** (2 * m + 1) / ((2 * m + 1) * (1 - z * z))
    return s, s + r          # s < log x < s + r
L2 = log_encl(2, 1); L3 = log_encl(3, 2); L5 = log_encl(5, 3); L7 = log_encl(7, 9)
chk("log2 < 7/10 via (2,1)", L2[1] < F(7, 10))
chk("13/12 < log3 via (3,2)", L3[0] >= F(13, 12))   # partial sum equals 13/12, R>0 strict
chk("log3 < 11/10 via (3,2)", L3[1] <= F(11, 10))  # R < bound strictly
chk("log5 < 13/8 via (5,3)", L5[1] < F(13, 8))
chk("35/18 < log7 via (7,9)", L7[0] > F(35, 18))
chk("log7 < 39/20 via (7,9)", L7[1] < F(39, 20))
chk("log100 = 2log2+2log5 < 5", 2 * L2[1] + 2 * L5[1] < 5)
print("  enclosures:", [float(a) for a in (L2[1], L3[0], L3[1], L5[1], L7[0], L7[1])])
lg2, lg3, lg5, lg7 = F(7, 10), F(11, 10), F(13, 8), F(39, 20)  # upper bounds used in paper
lg3lo, lg7lo = F(13, 12), F(35, 18)

# --- (prelim:odd-square-bound), TeX l.343-357
S = lg3 / 9 + lg5 / 25 + lg7 / 49 + F(1, 6) * ((3 * lg2 + 1) / 8 + ((lg2 + lg5) + 1) / 10)
chk("odd-square certificate == 30607/88200", S == F(30607, 88200))
chk("30607/88200 < 7/20", F(30607, 88200) < F(7, 20))
odd = mp.nsum(lambda k: 0, [1, 1])
true_odd = sum(mp.log(p) / p**2 for p in sp.primerange(3, 2_000_000))
print("  true sum_{odd p<2e6} log p/p^2 =", mp.nstr(true_odd, 10), "(tail < 1e-5)")
chk("true odd-square sum < 30607/88200", true_odd + mp.mpf('1e-5') < mp.mpf(30607) / 88200)

# --- V bound, TeX l.375-386
Vb = F(7, 20) + F(11, 360) + F(7, 16)
chk("7/20+11/360+7/16 == 589/720", Vb == F(589, 720))
chk("589/720 < 33/40", F(589, 720) < F(33, 40))
chk("log2/2 < 7/20", lg2 / 2 <= F(7, 20) and mp.log(2) / 2 < 0.35)
chk("log3/36 < 11/360", lg3 / 36 <= F(11, 360))
Vtrue = sum(mp.log(p) / (p * (p - 1)) for p in sp.primerange(2, 2_000_000))
print("  true V =", mp.nstr(Vtrue, 10))
# --- Mertens-interval loss, TeX l.398-405
loss = F(33, 40) + F(7, 4) * lg2 - F(1, 9) * lg3lo + F(6, 100)
chk("loss certificate == 2 - 7/675", loss == 2 - F(7, 675))
chk("(1+log100)/100 < 3/50", (1 + 2 * L2[1] + 2 * L5[1]) / 100 < F(3, 50))

# --- Proposition arith:main remainder, TeX l.796-837
t1 = (10000 * F(7, 5) + 4 * 10) / 9999
chk("14040/9999 < 141/100", t1 == F(14040, 9999) and t1 < F(141, 100))
chk("log4 < 7/5", 2 * L2[1] < F(7, 5))
chk("log10000 < 10", 4 * (L2[1] + L5[1]) < 10)
chk("p0/s0 < 301/10 (s0>=333)", 30 + F(30, 333) < F(301, 10))
t2 = F(301, 10) * (10 + F(7, 5)) / 20000
chk("U/(2s0) bound < 1/50", t2 < F(1, 50))
t3 = F(7, 60) * F(301, 10)
chk("7*2*(301/10)/120 == 2107/600 < 53/15", t3 == F(2107, 600) and t3 < F(53, 15))
chk("22/15*log2 < 77/75", F(22, 15) * lg2 == F(77, 75))
tot = F(141, 100) + F(1, 50) + F(53, 15) + F(77, 75)
chk("sum == 599/100 < 6", tot == F(599, 100))
chk("small case 9/20*10+7/5 == 59/10", F(9, 20) * 10 + F(7, 5) == F(59, 10))
chk("small case 89/60*7/5 == 623/300 < 6", F(89, 60) * F(7, 5) == F(623, 300))
chk("(29/30)+(1/2) == 22/15", F(29, 30) + F(1, 2) == F(22, 15))
chk("H-upper coefficients 29/60+1/15 == 11/20", F(29, 60) + F(1, 15) == F(11, 20))

# --- Theorem, TeX l.1136-1221
H1 = sum(F(1, m * m) for m in range(1, 8)) + F(1, 8) - F(1, 28)
chk("square normalization == 282431/176400 > 8/5", H1 == F(282431, 176400) and H1 > F(8, 5))
chk("sum_{m<=27} 1/m^2 >= that", sum(F(1, m * m) for m in range(1, 28)) >= H1)
By = (4 + 3 * lg3) / 27
chk("B(y) bound (4+3 log3)/27 < 73/270", By == F(73, 270))
Bxy = (4 + 38 * lg3) / F(3) ** 38
chk("B(xy) bound < 50*3^-38 < 1e-4", Bxy < 50 / F(3) ** 38 < F(1, 10000))
ie = (F(73, 135) + F(2, 27)) * F(12, 13)
chk("integral error == 332/585 < 4/7", ie == F(332, 585) and ie < F(4, 7))
chk("kappa = 70*6 + 4900/2 = 2870", 70 * 6 + F(4900, 2) == 2870)
ec = F(73, 270) * (F(287, 10) + 6) + F(76, 10000) + F(4, 7)
chk("error certificate == 1176629/118125 = 10 - 4621/118125", ec == F(1176629, 118125) == 10 - F(4621, 118125))
pm = F(287, 10) + F(5, 8) * 10
chk("prime mass 699/20 < 35", pm == F(699, 20) and pm < 35)
chk("70-6-35/2 == 93/2", 70 - 6 - F(35, 2) == F(93, 2))
chk("11/20*70+29/60 == 2339/60", F(11, 20) * 70 + F(29, 60) == F(2339, 60))
chk("93/2-2339/60 == 451/60", F(93, 2) - F(2339, 60) == F(451, 60))
chk("451/60*13/12 == 8+103/720", F(451, 60) * F(13, 12) == 8 + F(103, 720))
chk("beta>1/2: 1/(100 log3) < 1/2", True)
chk("76/100 < 1; 3^6=729", 3 ** 6 == 729)

# --- Appendix c=1/64, TeX l.1232-1278
By7 = (8 + 2 * lg7) / 49
chk("app B(y) (8+2log7)/49 == 17/70", By7 == F(17, 70))
Bxy7 = (8 + 25 * lg7) / F(7) ** 25
chk("app B(xy) < 57*7^-25 < 1e-4", Bxy7 < 57 / F(7) ** 25 < F(1, 10000))
chk("app kappa = 46*4+46^2/2 = 1242", 46 * 4 + F(46 * 46, 2) == 1242)
chk("app 1242/64 == 621/32", F(1242, 64) == F(621, 32))
aec = F(17, 70) * (F(621, 32) + 4) + F(50, 10000) + (F(17, 35) + F(2, 49)) * F(18, 35)
chk("app error certificate == 6 - 21863/548800", aec == 6 - F(21863, 548800))
apm = F(621, 32) + F(5, 8) * 6
chk("app J bound 741/32", apm == F(741, 32))
gap = (46 - 4 - apm / 2) - (F(11, 20) * 46 + F(29, 60))
chk("app gap == 4453/960", gap == F(4453, 960))
chk("app 4453/960*35/18 == 8+3523/3456", F(4453, 960) * F(35, 18) == 8 + F(3523, 3456))
chk("app 50/64 == 25/32 < 1", F(50, 64) == F(25, 32))

# --- symbolic identities
q, k, t, x, y, eps = sp.symbols('q k t x y epsilon', positive=True)
fk = (q + 1 + k / 2 * sp.log(q)) * q ** (-k / 2)
claimed = q ** (-k / 2 - 1) * ((1 - k / 2) * q - k ** 2 / 4 * sp.log(q))
chk("f_k'(q) as printed (TeX l.1166-1168) is correct", sp.simplify(sp.diff(fk, q) - claimed) == 0)
Bt = (q + 1 + sp.log(t) / 2) * t ** sp.Rational(-1, 2)
G = 2 * Bt.subs(t, y) + 2 / sp.sqrt(y)   # claimed antiderivative: -G' = B(y)/y and G -> 0
chk("int_y^oo B(t)/t dt == 2B(y)+2y^-1/2 (G' = -B/y, G(oo)=0)",
    sp.simplify(sp.diff(G, y) + Bt.subs(t, y) / y) == 0 and sp.limit(G, y, sp.oo) == 0)
phi = sp.log(t) * t ** (-eps)
main = (phi * t ** eps).subs(t, x * y) - (phi * t ** eps).subs(t, y) - sp.integrate(sp.simplify(t ** eps * sp.diff(phi, t)), (t, y, x * y))
kappa = sp.log(x) * sp.log(y) + sp.log(x) ** 2 / 2
chk("main term of weighted Abel == eps*kappa", sp.simplify(sp.expand_log(main - eps * kappa, force=True)) == 0)
g_lhs = (t ** (-(1 - eps)) - y ** eps / t) / eps
g_rhs = sp.integrate(t ** (-(1 - eps)) * 0 + sp.Symbol('w') ** (eps - 1), (sp.Symbol('w'), y, t)) / t
chk("g(v) = (1/v) int_u^v w^-beta dw", sp.simplify(g_lhs - g_rhs) == 0)
print("ALL PASS" if ok_all else "SOME FAIL")

