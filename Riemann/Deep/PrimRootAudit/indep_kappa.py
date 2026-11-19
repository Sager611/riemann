#!/usr/bin/env python3
"""indep_kappa.py (ppraudit) -- independent recomputation of Theorem A entries (DERIVATION sec 3-4) for selected omega.
Own code: explicit-stack enumeration of Rosser's set D(y) (no bulk leaves; T in units 2^-128), kernel constants and Lemma Z
integrals by mpmath quadrature (not EP's closed forms), sup_{L>=Lm} Z(L9)/L by a dense scan, margin (4.4) at 256 bits.
Usage: python3 indep_kappa.py w1 w2 ...   (levels y_w read from PrimRootNum/logs/sieve_best.json; any y is valid)"""
import sys, json, math, time
from fractions import Fraction as Fr
import mpmath as M
M.mp.prec = 256
B = M.mpf(8733819666132076) / 10 ** 16; b = 1 - B; IB = 1 / b
ell, cpk = M.mpf(7), M.mpf(21) / 10; a_ = ell - cpk           # rise 4.9, fall 2.1
def g(u): return (u + ell) / a_ if u <= -cpk else -u / cpk
def I(fun): return M.quad(fun, [-ell, -cpk, 0])
W1 = I(lambda u: g(u) * M.e ** (b * u)); W0 = I(lambda u: g(u) * M.e ** (-B * u))
winf = max(M.findroot(lambda u: 0, 0) * 0 + max(g(u) * M.e ** (-B * u) for u in [-ell + (ell) * k / 20000 for k in range(20001)]),
           ((1 / B) / a_) * M.e ** (B * (ell - 1 / B)))          # exact max of rising part at u = -7 + 1/B
m = M.mpf(1) / 2 - B
n0, n1, n2 = ell / 2, M.mpf(2), 2 * (1 / a_ + 1 / cpk)
m0 = I(lambda u: g(u) * M.e ** (m * u))
J1 = I(lambda u: (1 / a_ if u < -cpk else 1 / cpk) * M.e ** (m * u))
J2 = M.e ** (-m * ell) / a_ + (1 / a_ + 1 / cpk) * M.e ** (-m * cpk) + 1 / cpk
m1 = J1 + abs(m) * m0; m2 = J2 + 2 * abs(m) * J1 + m * m * m0
SUMB = 1 + M.euler / 2 - M.log(4 * M.pi) / 2; SIG = 2 * (1 + B * B / 196) * SUMB
CB = M.mpf(22737) / 100000; S0 = M.mpf(5) / 7; UE = 2 * M.pi - 2
aa = lambda t: t / M.pi + CB
bb = lambda t: (t / M.pi) * M.log(t / (2 * M.pi * M.e)) + M.mpf(1) / 4 + CB * M.log((t + 2) / (2 * M.pi)) - M.mpf(1) / 2
rho0 = lambda t: M.log(t / (2 * M.pi)) / M.pi + CB / (t + 2)
def zco(N0, K1, K2):
    G = lambda t: min(N0, K1 / t, K2 / t ** 2)
    br = sorted(set([S0, max(S0, K1 / N0), max(S0, K2 / K1), max(S0, M.sqrt(K2 / N0)), max(S0, UE)]))
    pts = br + [M.inf]
    za = aa(S0) * N0 + M.quad(lambda t: G(t), pts) / M.pi
    zb = bb(S0) * N0 + M.quad(lambda t: rho0(t) * G(t), pts)
    pe = [p for p in pts if p >= max(S0, UE)]
    ze = 2 * M.quad(lambda t: G(t) / (t + 2), pe)
    return za, zb, ze
def Z(co, N0, L): za, zb, ze = co; return L * za + 2 * N0 * M.log(1 + L) + zb + ze / (1 + L)
L9 = M.log(9)
def supZL(co, N0, Lm):
    best = M.mpf(0); L = Lm
    for k in range(1600):           # L = Lm e^{k/40}: up to Lm e^40; beyond, Z/L -> za (kept by the final max)
        best = max(best, Z(co, N0, max(L, L9)) / L); L = L * M.e ** (M.mpf(1) / 40)
    return max(best, co[0])
def margin(kap, Lm, T, C):
    X = M.log(kap * Lm) * IB; Y = X - ell
    if Y <= 1: return None
    eY = M.e * Y
    k1, k2 = n1 + n0 / eY, n2 + 2 * n1 / eY + (2 / eY) ** 2 * n0
    k1L, k2L = m1 + m0 / eY, m2 + 2 * m1 / eY + (2 / eY) ** 2 * m0
    ZoL = supZL(zco(n0, k1, k2), n0, Lm); ZLoL = supZL(zco(m0, k1L, k2L), m0, Lm)
    xB, xh = M.e ** (B * X), M.e ** ((M.mpf(1) / 2 - B) * X)
    z1, z1L = SIG * k2 / 2, SIG * k2L / 2
    triv = W0 / ((1 - M.e ** (-2 * Y)) * xB * Lm)
    Xo = W1 * kap - (z1 + xh * z1L) / Lm - triv - winf * (7 / Lm + 1) / xB
    Ro = ZoL / 2 + xh * ZLoL / 2 + triv
    return T * Xo - (C - 1) * Ro
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]
PR = primes(120000)
K = 128; ONE = 1 << K
def rosser(w, y):
    """explicit DFS over D(y) restricted to d | P_w; returns (T lower bound numerator 2^-128, C, exact T if small)"""
    P = PR[:w]; Tn = ONE; C = 1; exact = Fr(1) if w <= 14 else None
    st = [(1, w, 0)]                  # (prod, number of admissible indices j < li, length r)
    while st:
        prod, li, r = st.pop()
        mpos = r + 1
        for j in range(li):
            q = P[j]
            if mpos % 2 == 0 and prod * q * q * q >= y: break
            d = prod * q
            if mpos % 2 == 1: Tn -= -((-ONE) // d)          # mu = -1: subtract ceil
            else: Tn += ONE // d                             # mu = +1: add floor
            C += 1
            if exact is not None: exact += Fr((-1) ** mpos, d)
            if j > 0: st.append((d, j, mpos))
    return Tn, C, exact
tab = {}
for line in open("../PrimRootNum/logs/kappa_table.txt"):
    if line.startswith("#"): continue
    f = line.split(); tab[int(f[0])] = (int(f[2]), float(f[3]), int(f[4]), float(f[5]), Fr(f[6]), float(f[7]))
print("kernel: W1 = %s W0 = %s winf = %s m0 = %s m1 = %s m2 = %s SIG = %s" % tuple(M.nstr(v, 10) for v in (W1, W0, winf, m0, m1, m2, SIG)), flush=True)
for w in [int(a) for a in sys.argv[1:]]:
    t0 = time.time(); y, Ttab, Ctab, Lmtab, kap, cw = tab[w]
    Tn, C, ex = rosser(w, y)
    T = M.mpf(Tn) / ONE
    Pw = 1
    for q in PR[:w]: Pw *= q
    Lm = M.log(Pw + 1) if Pw + 1 >= 3 else M.log(3)
    mg = margin(M.mpf(kap.numerator) / kap.denominator, Lm, T, C)
    # least kappa with positive margin (own), 3 significant digits; bisection at 64 bits, endpoint re-checked at 256 bits
    M.mp.prec = 64
    lo_, hi_ = M.mpf(kap.numerator) / kap.denominator / 50, M.mpf(kap.numerator) / kap.denominator * 2
    for _ in range(24):
        mid = M.sqrt(lo_ * hi_); mm = margin(mid, Lm, T, C)
        if mm is not None and mm > 0: hi_ = mid
        else: lo_ = mid
    M.mp.prec = 256
    e = int(M.floor(M.log10(hi_))) - 2; kmin = Fr(int(M.ceil(hi_ / M.mpf(10) ** e)), 1) * Fr(10) ** e
    c3 = M.cbrt(M.mpf(kap.numerator) / kap.denominator / Lm ** 2)
    print("w=%d y=%d: C %d (table %d) %s; T %.10g (table %.8g)%s; Lm %.6f (table %.6f); kappa_table %s margin %s; own least kappa %s (%s); c_w %.6f (table %.6f) [%.0f s]" % (
        w, y, C, Ctab, "OK" if C == Ctab else "MISMATCH", float(T), Ttab, (" exact %.12g" % float(ex)) if ex is not None else "",
        float(Lm), Lmtab, kap, M.nstr(mg, 8), kmin, "<= table OK" if kmin <= kap else "ABOVE TABLE", float(c3), cw, time.time() - t0), flush=True)
