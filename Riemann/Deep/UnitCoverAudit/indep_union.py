#!/usr/bin/env python3
"""indep_union.py (unitaudit) -- independent 256-bit INTERVAL (mpmath.iv) recomputation of UnitCover sec 3.4-5:
a1, b1, the certified sup over Q of the union bound U (linear theorem), the table rows U(Q) < 1 at the printed (r, s, alpha),
local minimality of the rows, and U(Q) at the actual ceilings on a grid of Q.  Classes from logs/classes_k*.json
(own enumeration, indep_classes.py)."""
import json, os, math
from fractions import Fraction as Fr
from mpmath import iv, mp
iv.prec = 256; mp.prec = 256
HERE = os.path.dirname(os.path.abspath(__file__))
CL = {k: [(Fr(p), m, n, D, c) for (p, m, n, D, c) in json.load(open(os.path.join(HERE, "logs", "classes_k%d.json" % k)))] for k in (1, 2, 3, 4)}
def I(x):
    x = Fr(x); return iv.mpf(x.numerator) / iv.mpf(x.denominator)
def hi(x): return mp.mpf(x.b)
def lo(x): return mp.mpf(x.a)
def params(k, d, al):
    d, al = I(d), I(al)
    beta = d + (1 - d) * al; t = al / (1 - al)
    rg = I(Fr(1, 2 ** k)) * (1 - d) + d
    return d, al, beta, t, (k + 1) / iv.log(1 / rg), 2 / iv.log(1 / beta)
def supf(f, eps, Lmin):
    if f == 0: return iv.exp(-eps * Lmin)
    Ls = f / eps - 1
    atmin = (1 + Lmin) ** f * iv.exp(-eps * Lmin)
    glob = (f / eps) ** f * iv.exp(eps - f)
    if hi(Ls) <= lo(Lmin): return atmin
    return iv.mpf([0, max(hi(atmin), hi(glob))]) if lo(Ls) < hi(Lmin) else glob
def linear(k, d, al, a0, b0, Lmin):
    d, al, beta, t, a1, b1 = params(k, d, al)
    tot = iv.mpf(0)
    for idx, (p, m, n, D, c) in enumerate(CL[k]):
        rho = I(p) + d * (1 - I(p))
        e, f, coef = ((1, 0, iv.mpf(1)) if m == 0 else (m + 1, 0, 1 / iv.mpf(m + 1))) if n == 0 else (m + 1, D - 1, iv.mpf(2) ** n)
        if idx == 0:
            assert (p, m, n, D, c) == (Fr(1, 2 ** k), k, 0, 1, 1); mfac = iv.mpf(1)
        else:
            eps = a1 * iv.log(1 / rho) - e
            assert lo(eps) > 0, (k, p, m)
            mfac = supf(f, eps, Lmin)
        tot += c * coef * mfac * rho ** a0
    return a1, b1, tot / t ** k + al / 2 * beta ** b0
def U(Q, k, d, al, r, s):
    d, al, beta, t, a1, b1 = params(k, d, al)
    tot = iv.mpf(0)
    for (p, m, n, D, c) in CL[k]:
        rho = I(p) + d * (1 - I(p))
        qs = (Q if m == 0 else Q ** (m + 1) / (m + 1)) if n == 0 else iv.mpf(2) ** n * Q ** (m + 1) * (1 + iv.log(Q)) ** (D - 1)
        tot += c * qs * rho ** r
    return tot / t ** k + al * Q ** 2 / 2 * beta ** s
L3 = iv.log(3)
print("== linear theorem (page / D:136-147 values in brackets)")
for (k, d, al, a0, Lmin, claim) in [(1, Fr(1, 1000), Fr(1, 10**4), 13, L3, "0.6183"), (2, Fr(1, 1000), Fr(1, 10**4), 29, L3, "0.4641"),
                                    (3, Fr(1, 1000), Fr(1, 10**4), 42, L3, "0.7266"), (4, Fr(1, 1000), Fr(1, 10**4), 56, L3, "0.5464"),
                                    (4, Fr(1, 1000), Fr(1, 10**4), 13, 100 * iv.log(10), "0.5388"),
                                    (3, Fr(1, 10**6), Fr(1, 10**7), 72, L3, "0.6524"), (4, Fr(1, 10**6), Fr(1, 10**7), 96, L3, "0.4770")]:
    a1, b1, tot = linear(k, d, al, a0, 0, Lmin)
    a1m, b1m, totm = linear(k, d, al, a0 - 1, 0, Lmin)
    print("k=%d delta=%s alpha=%s a0=%d: a1 <= %s, b1 <= %s, a1+b1 <= %s, sup U <= %s [%s]; with a0-1: %s"
          % (k, d, al, a0, mp.nstr(hi(a1), 8), mp.nstr(hi(b1), 8), mp.nstr(hi(a1 + b1), 8), mp.nstr(hi(tot), 6), claim, mp.nstr(hi(totm), 4)), flush=True)
print("1/log 2 + 0.664 = %s; 5/(4 log 2) = %s; 2.1067*log 2 = %s" % (mp.nstr(1 / mp.log(2) + mp.mpf('0.664'), 8), mp.nstr(5 / (4 * mp.log(2)), 6), mp.nstr(mp.mpf('2.1067') * mp.log(2), 6)))
print("== table rows: U(Q) at the printed (k, r, s, alpha); also (r-1, s), (r, s-1)")
rows = [(3, Fr(1, 1000), [(1, 25, 3, Fr(3, 100)), (2, 21, 6, Fr(1, 10)), (3, 24, 6, Fr(1, 10)), (4, 27, 6, Fr(1, 10))]),
        (6, Fr(1, 1000), [(1, 46, 6, Fr(1, 100)), (2, 37, 7, Fr(3, 100)), (3, 35, 12, Fr(1, 10)), (4, 38, 12, Fr(1, 10))]),
        (12, Fr(1, 1000), [(1, 88, 9, Fr(3, 1000)), (2, 67, 12, Fr(1, 100)), (3, 62, 15, Fr(3, 100)), (4, 67, 15, Fr(3, 100))]),
        (30, Fr(1, 1000), [(1, 209, 22, Fr(1, 1000)), (2, 161, 20, Fr(3, 10**4)), (3, 143, 22, Fr(1, 1000)), (4, 134, 30, Fr(1, 100))]),
        (100, Fr(1, 1000), [(1, 678, 67, Fr(1, 10**4)), (2, 512, 67, Fr(1, 10**4)), (3, 458, 67, Fr(1, 10**4)), (4, 431, 67, Fr(1, 10**4))]),
        (1000, Fr(1, 1000), [(1, 6670, 666, Fr(1, 10**5)), (2, 5010, 666, Fr(1, 10**5)), (3, 4461, 666, Fr(1, 10**5)), (4, 4192, 666, Fr(1, 10**5))]),
        (12, Fr(1, 10**6), [(3, 62, 15, Fr(3, 100)), (4, 67, 15, Fr(3, 100))]),
        (100, Fr(1, 10**6), [(3, 464, 33, Fr(3, 10**7)), (4, 435, 35, Fr(1, 10**6))]),
        (1000, Fr(1, 10**6), [(3, 4452, 335, Fr(1, 10**7)), (4, 4177, 333, Fr(3, 10**8))])]
for (e10, d, lst) in rows:
    Q = iv.mpf(10) ** e10; out = []
    for (k, r, s, al) in lst:
        u = U(Q, k, d, al, r, s); u1 = U(Q, k, d, al, r - 1, s); u2 = U(Q, k, d, al, r, s - 1) if s > 0 else iv.mpf(9)
        ok = hi(u) < 1; loc = lo(u1) >= 1 and lo(u2) >= 1
        out.append("k=%d N=%d U<=%s %s%s" % (k, r + s, mp.nstr(hi(u), 4), "ok" if ok else "FAIL", "" if loc else " (not locally minimal: %s %s)" % (mp.nstr(lo(u1), 4), mp.nstr(lo(u2), 4))))
    print("Q=10^%d delta=%s: %s" % (e10, d, "; ".join(out)), flush=True)
print("== actual U(Q) at r = ceil(a1 log Q + 56), s = ceil(b1 log Q), k=4, delta=1e-3, alpha=1e-4, on a grid")
d, al, beta, t, a1, b1 = params(4, Fr(1, 1000), Fr(1, 10**4))
worst = (0, None)
for lq in [math.log(3), math.log(4), math.log(9), 3, 5, 10, 20, 50, 100, 230, 500, 1000, 2303, 10**4, 10**5]:
    Q = iv.exp(iv.mpf(lq)); r = int(math.ceil(float(mp.mpf((a1 * lq + 56).b)))); s = int(math.ceil(float(mp.mpf((b1 * lq).b))))
    u = U(Q, 4, Fr(1, 1000), Fr(1, 10**4), r, s)
    if hi(u) > worst[0]: worst = (hi(u), lq)
    print("  log Q = %g: r = %d, s = %d, N = %d (2.1067 log Q + 58 = %.2f), U <= %s" % (lq, r, s, r + s, 2.1067 * lq + 58, mp.nstr(hi(u), 4)))
print("max on grid %s at log Q = %g" % (mp.nstr(worst[0], 4), worst[1]))
