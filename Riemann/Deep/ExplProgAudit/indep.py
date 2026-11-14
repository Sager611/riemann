#!/usr/bin/env python3
"""indep.py -- independent recomputation (explaudit) of the ExplProg constants at 256 bits with mpmath.
Own implementation from DERIVATION.md formulas; zero sums by numerical quadrature (mp.quad), NOT the closed-form
antiderivatives of cert.py; kernel constants by quadrature. Not rigorous interval arithmetic: margins are compared
with cert.py's certified values. Modes: const | lemz | linnik | nonres | pnt <cells-file> | ptail | lprime | small"""
import sys, time, math
import mpmath as mp
mp.mp.prec = int(__import__("os").environ.get("INDEP_PREC", "256"))
pi = mp.pi; E = mp.e; EUL = mp.euler
B2 = (777 - mp.sqrt(1081)) / 852
B = mp.mpf(8733819666132076) / 10**16
b = 1 - B; INVB = 1 / b
CB = mp.mpf(22737) / 100000
S0 = mp.mpf(5) / 7
HALF = mp.mpf(1) / 2
LOG2 = mp.log(2); LOG9 = mp.log(9)
UE = 2 * pi - 2

def env(n0, k1, k2):
    if k1 * k1 <= n0 * k2:
        ta = k1 / n0; tb = k2 / k1
        return [(mp.mpf(0), ta, n0, 0), (ta, tb, k1, 1), (tb, mp.inf, k2, 2)]
    tc = mp.sqrt(k2 / n0)
    return [(mp.mpf(0), tc, n0, 0), (tc, mp.inf, k2, 2)]

def qpieces(pcs, f, lo):
    s = mp.mpf(0)
    for (a, bb, c, m) in pcs:
        a2 = a if a > lo else lo
        if bb <= a2: continue
        s += c * mp.quad(lambda t: f(t) * t ** (-m), [a2, bb])
    return s

def rho0(t): return mp.log(t / (2 * pi)) / pi + CB / (t + 2)
A57 = S0 / pi + CB
B57 = (S0 / pi) * mp.log(S0 / (2 * pi * E)) + mp.mpf(1) / 4 + CB * mp.log((S0 + 2) / (2 * pi)) - HALF

def zco(n0, k1, k2):
    pcs = env(n0, k1, k2)
    za = A57 * n0 + qpieces(pcs, lambda t: 1 / pi, S0)
    zb = B57 * n0 + qpieces(pcs, rho0, S0)
    ze = 2 * qpieces(pcs, lambda t: 1 / (t + 2), UE if UE > S0 else S0)
    return za, zb, ze

def Zof(co, n0, L):
    za, zb, ze = co
    return L * za + 2 * n0 * mp.log(1 + L) + zb + ze / (1 + L)

def Rbmor(t, L):  # undecomposed BMOR upper bound (chi(-1) term as +1/4), conductor e^L
    l = L + mp.log((t + 2) / (2 * pi))
    return (t / pi) * (L + mp.log(t / (2 * pi * E))) + mp.mpf(1) / 4 + CB * l + 2 * mp.log(1 + l) - HALF

def stieltjes_bmor(n0, k1, k2, L):
    """sum_rho G <= R(5/7) n0 + int_{5/7}^oo R'(t) G(t) dt with the undecomposed R (independent of the N+ split)"""
    pcs = env(n0, k1, k2)
    dR = lambda t: mp.diff(lambda s: Rbmor(s, L), t)
    return Rbmor(S0, L) * n0 + qpieces(pcs, dR, S0)

CR0, CR1, CR2 = mp.mpf(4926) / 1000, mp.mpf(274) / 1000, mp.mpf(886) / 1000
def Nz(t): return (t / pi) * mp.log(t / (2 * pi * E)) + CR0 + CR1 * mp.log(t) + CR2 * mp.log(mp.log(t))
def zeta_sum(n0, k1, k2):
    pcs = env(n0, k1, k2); T = mp.mpf(14)
    G14 = min(n0, k1 / T, k2 / T ** 2)
    rz = lambda t: mp.log(t / (2 * pi)) / pi + CR1 / t + CR2 / (t * mp.log(t))
    return Nz(T) * G14 + qpieces(pcs, rz, T)

SUM_BETA = 1 + EUL / 2 - mp.log(4 * pi) / 2
SIG = 2 * (1 + B * B / 196) * SUM_BETA

# ---- Linnik kernel by quadrature
ELL = mp.mpf(7); RISE = mp.mpf(49) / 10; CC = ELL - RISE
def g(u): return (u + ELL) / RISE if u <= -CC else -u / CC
def gp(u): return 1 / RISE if u < -CC else -1 / CC
def kq(f): return mp.quad(f, [-ELL, -CC, 0])
W1 = kq(lambda u: g(u) * mp.exp(b * u))
W0 = kq(lambda u: g(u) * mp.exp(-B * u))
ustar = -ELL + 1 / B
WINF = max(g(ustar) * mp.exp(-B * ustar), g(-CC) * mp.exp(B * CC))
n0, n1, n2 = ELL / 2, mp.mpf(2), 2 * (1 / RISE + 1 / CC)
mL = HALF - B
m0 = kq(lambda u: g(u) * mp.exp(mL * u))
J1 = kq(lambda u: abs(gp(u)) * mp.exp(mL * u))
J2 = mp.exp(-mL * ELL) / RISE + (1 / RISE + 1 / CC) * mp.exp(mL * (-CC)) + 1 / CC
m1 = J1 + abs(mL) * m0
m2 = J2 + 2 * abs(mL) * J1 + mL * mL * m0
m1true = kq(lambda u: abs(gp(u) * mp.exp(mL * u) + mL * g(u) * mp.exp(mL * u)))
def ks(Y):
    eY = E * Y
    return (n0, n1 + n0 / eY, n2 + 2 * n1 / eY + (2 / eY) ** 2 * n0), (m0, m1 + m0 / eY, m2 + 2 * m1 / eY + (2 / eY) ** 2 * m0)

def s(x, d=10): return mp.nstr(x, d)

def dRbmor(t, L):
    l = L + mp.log((t + 2) / (2 * pi))
    return (L + mp.log(t / (2 * pi))) / pi + CB / (t + 2) + 2 / ((1 + l) * (t + 2))

def stieltjes_bmor2(n0, k1, k2, L):
    pcs = env(n0, k1, k2)
    return Rbmor(S0, L) * n0 + qpieces(pcs, lambda t: dRbmor(t, L), S0)

def zmaj(co, nn, Lm):
    za, zb, ze = co
    cond = 2 * nn * (mp.log(1 + Lm) - Lm / (1 + Lm)) + zb
    zbb = zb if cond >= 0 else max(zb, 0)
    return za + (2 * nn * mp.log(1 + Lm) + zbb) / Lm + max(ze, 0) / (Lm * (1 + Lm))

def F_linnik(c2, Lm, Ldir=None, phi=False):
    """tail majorant at Lm (Ldir None) or direct value of (4.1)/(q x^B L) at L = Ldir with k's at the true Y"""
    L = Lm if Ldir is None else Ldir
    qv = mp.exp(L)
    if phi:
        qv = qv / (mp.exp(EUL) * mp.log(L) + mp.mpf(250637) / 100000 / mp.log(L))
    X = mp.log(c2 * qv * L) / b; Y = X - ELL
    (a0, k1, k2), (c0, k1L, k2L) = ks(Y)
    co = zco(a0, k1, k2); coL = zco(c0, k1L, k2L)
    xhB = mp.exp((HALF - B) * X); xB = mp.exp(B * X); x3B = mp.exp((mp.mpf(1) / 3 - B) * X)
    if Ldir is None:
        main = zmaj(co, a0, L) / 2 + xhB * zmaj(coL, c0, L) / 2
    else:
        Lq = max(L, LOG9)
        main = (Zof(co, a0, Lq) + xhB * Zof(coL, c0, Lq)) / (2 * L)
    rest = SIG * (k2 + xhB * k2L) / 2 / (qv * L) + W0 / ((1 - mp.exp(-2 * Y)) * xB * L) \
        + WINF * (ELL / LOG2 + 1) / xB + mp.mpf(101624) / 100000 * WINF * (xhB + x3B * X / LOG2) / L
    lhs = c2 * W1
    return (lhs - main - rest) / lhs

def F_nonres(c3, Lm, Ldir=None):
    L = Lm if Ldir is None else Ldir
    X = mp.log(c3 * L) / b; Y = X - ELL
    (a0, k1, k2), (c0, k1L, k2L) = ks(Y)
    co = zco(a0, k1, k2); coL = zco(c0, k1L, k2L)
    xhB = mp.exp((HALF - B) * X); xB = mp.exp(B * X)
    if Ldir is None:
        main = zmaj(co, a0, L) / 2 + xhB * zmaj(coL, c0, L) / 2
    else:
        Lq = max(L, LOG9)
        main = (Zof(co, a0, Lq) + xhB * Zof(coL, c0, Lq)) / (2 * L)
    rest = SIG * (k2 + xhB * k2L) / 2 / L + 2 * W0 / ((1 - mp.exp(-2 * Y)) * xB * L) + 2 * WINF * (ELL / LOG2 + 1) / xB
    lhs = c3 * W1
    return (lhs - main - rest) / lhs

def mode_const():
    print("B2 =", s(B2, 22), " B - B2 =", s(B - B2, 5), " (need > 2.2e-254)", "  1/b =", s(INVB, 14))
    print("a(5/7) =", s(A57), " b(5/7) =", s(B57), "  R(5/7;9) =", s(Rbmor(S0, LOG9)), "  SIG =", s(SIG, 12), " SIG/2 =", s(SIG / 2, 12))
    print("W(1) =", s(W1), " W(0) =", s(W0), " ||w||inf =", s(WINF), " n2 =", s(n2))
    print("m0 =", s(m0), " m1 =", s(m1), " (true ||gL'||_1 =", s(m1true), ") m2 =", s(m2))
    co = zco(n0, n1, n2)
    print("Y=oo: za =", s(co[0]), " c_inf = za/(2W1) =", s(co[0] / (2 * W1)))
    X = INVB * mp.log(10**4 * mp.log(10**4)); Y = X - ELL
    (a0, k1, k2), _ = ks(Y)
    co = zco(a0, k1, k2)
    print("q=1e4, x=(q log q)^(1/b): Y =", s(Y, 6), " k1 =", s(k1, 7), " k2 =", s(k2, 7), " za,zb,ze =", [s(v, 7) for v in co])
    print("share of za from zeros below 5/7:", s(A57 * a0 / co[0], 4))

def mode_lemz():
    for L in [LOG9, mp.log(10**4), mp.mpf(50)]:
        X = INVB * mp.log(mp.exp(L) * L); Y = X - ELL
        (a0, k1, k2), (c0, k1L, k2L) = ks(Y)
        Zc = Zof(zco(a0, k1, k2), a0, L); Zs = stieltjes_bmor2(a0, k1, k2, L)
        ZcL = Zof(zco(c0, k1L, k2L), c0, L); ZsL = stieltjes_bmor2(c0, k1L, k2L, L)
        print("L=%s: Z(decomposed) = %s >= Stieltjes(undecomposed BMOR) = %s : %s;  left: %s >= %s : %s"
              % (s(L, 6), s(Zc), s(Zs), Zc >= Zs, s(ZcL), s(ZsL), ZcL >= ZsL))
    t = mp.mpf(3); L = mp.mpf(5)
    print("dR check:", s(dRbmor(t, L) - mp.diff(lambda u: Rbmor(u, L), t), 5))

def mode_linnik():
    rows = [(1, 2), ("0.79", 2), ("0.75", 3), ("0.72", 4), ("0.67", 6), ("0.62", 10), ("0.56", 20), ("0.53", 30), ("0.51", 50), ("0.49", 100), ("0.46", 1000)]
    for (c, e) in rows:
        c2 = mp.mpf(c); Lm = mp.floor(mp.log(mp.mpf(10) ** e) * 10**6) / 10**6
        print("Thm2 q>=1e%d c2=%s: tail margin %s (cert %s); phi-version %s" % (e, c, s(F_linnik(c2, Lm), 5), "see log", s(F_linnik(c2, Lm, phi=True), 5)), flush=True)
    # direct (non-frozen, non-majorised) values on a grid, c2 = 1 and 0.72
    worst = {}
    for c, Lstart in [(1, mp.log(100)), ("0.79", mp.log(100)), ("0.72", mp.log(10**4))]:
        c2 = mp.mpf(c); w = None
        for L in [Lstart * (mp.mpf(1.05) ** k) for k in range(0, 60, 3)]:
            v = F_linnik(c2, None, Ldir=L)
            w = v if w is None or v < w else w
        print("Thm2 direct grid c2=%s from L=%s: min margin %s" % (c, s(Lstart, 6), s(w, 5)), flush=True)

def mode_nonres():
    for (c, Ls) in [(1, 25), ("0.62", 25), ("0.60", 30), ("0.56", 50), ("0.52", 100), ("0.47", 1000), ("0.46", 10000)]:
        print("Thm3 log q>=%s c3=%s: tail margin %s" % (Ls, c, s(F_nonres(mp.mpf(c), mp.mpf(Ls)), 5)), flush=True)
    w = None
    for L in [mp.mpf(25) * (mp.mpf(1.07) ** k) for k in range(0, 60, 4)]:
        v = F_nonres(mp.mpf(1), None, Ldir=L); w = v if w is None or v < w else w
    print("Thm3 direct grid c3=1, L in [25, 25*1.07^56]: min margin", s(w, 5))
    # trivial range: q - 1 <= (log q)^(1/b) for integers 3 <= q <= e^25.6: concave check + endpoints + dense integers
    worst = None
    for q in list(range(3, 200000)):
        d = INVB * mp.log(mp.log(q)) - mp.log(q - 1)
        worst = d if worst is None or d < worst else worst
    print("Thm3 trivial: min over integers 3..2e5 of (1/b)loglog q - log(q-1) =", s(worst, 5), "; at log q=25.6:", s(INVB * mp.log(mp.mpf('25.6')) - mp.log(mp.exp(mp.mpf('25.6')) - 1), 5))

def mode_small():
    N = 70
    lam = [mp.mpf(0)] * (N + 1)
    for n in range(2, N + 1):
        p = next(d for d in range(2, n + 1) if n % d == 0); m = n
        while m % p == 0: m //= p
        if m == 1: lam[n] = mp.log(p)
    best = (0, None); best64 = (0, None)
    for q in range(1, 66):
        ph = sum(1 for k in range(1, q + 1) if math.gcd(k, q) == 1)
        for a in range(1, q + 1):
            if math.gcd(a, q) != 1: continue
            psi = mp.mpf(0)
            for n in range(1, 66):
                if n % q == a % q: psi += lam[n]
                if n < max(2, q): continue
                K = 400 if n < 8 else 8
                for j in range(K + 1):
                    x = n + mp.mpf(j) / K if j < K else n + 1 - mp.mpf(10) ** -30
                    r = abs(psi - x / ph) / (x ** B * mp.log(q * x) ** 2)
                    if r > best[0]: best = (r, (q, a, float(x)))
                    if x >= 64 and r > best64[0]: best64 = (r, (q, a, float(x)))
    print("exact small range sampled: max R =", s(best[0], 8), best[1], "; on [64,66):", s(best64[0], 6), best64[1])
    print("q=1, x=2: R =", s((2 - mp.log(2)) / (2 ** B * mp.log(2) ** 2), 8))

PSI_RS = mp.mpf(103883) / 100000
def trap(d, al):
    WB = ((1 + d) ** (B + 1) - 1) / (d * B * (B + 1)); A1 = (1 + d) ** B * (1 + al ** B); A2 = ((1 + d) ** (B + 1) + 1) * (1 + al ** B) / d
    Wh = ((1 + d) ** mp.mpf(1.5) - 1) / (d * mp.mpf(0.75)); A1L = (1 + d) ** HALF * (1 + al ** HALF); A2L = ((1 + d) ** mp.mpf(1.5) + 1) * (1 + al ** HALF) / d
    return (WB, A1, A2), (Wh, A1L, A2L)

def phil(L0, L1):
    return max(mp.exp(L0 / 2) / mp.sqrt(2), mp.exp(L0) / (mp.exp(EUL) * mp.log(L1) + mp.mpf(250637) / 100000 / mp.log(L0)))

def my_cell(kind, L0, L1, X0, X1, thetas=(0.25, 1, 4, 16, 64)):
    if kind == "gen": ph = phil(L0, L1); den = (L0 + X0) ** 2
    else: ph = mp.mpf(1); den = (X0 if kind == "q1" else LOG2 + X0) ** 2
    al = mp.exp(-X0 / 2); xh = mp.exp((HALF - B) * X0); best = None
    for th in thetas:
        d = min(HALF, th * ph * mp.exp(-b * X1))
        G, GL = trap(d, al)
        T1 = PSI_RS * xh
        T2 = d * (mp.exp(b * X1) + xh) / (2 * ph)
        T3 = (zeta_sum(*G) + xh * zeta_sum(*GL)) / (2 * ph)
        lg = mp.log((1 + d) / (1 - d))
        T5 = (X1 / 2 + lg) / (1 - (mp.exp(X0 / 2) * (1 - d)) ** -2) / mp.exp(B * X0)
        if kind == "gen":
            Lq = max(L1, LOG9)
            T4 = (Zof(zco(*G), G[0], Lq) + xh * Zof(zco(*GL), GL[0], Lq)) / 2
            T6 = ((L1 / LOG2) * (X1 / 2 + lg) + L1) / mp.exp(B * X0)
        elif kind == "q2": T4 = 0; T6 = (X1 / 2 + lg + LOG2) / mp.exp(B * X0)
        else: T4 = 0; T6 = 0
        R = (T1 + T2 + T3 + T4 + T5 + T6) / den
        best = R if best is None or R < best else best
    return best

def mode_pnt():
    import random
    rows = [l.split() for l in open(sys.argv[2])]
    random.seed(20261114)
    pick = sorted(range(len(rows)), key=lambda i: -float(rows[i][5]))[:3]
    pick += [i for i, r in enumerate(rows) if r[0] in ("q1", "q2")][:2] + [i for i, r in enumerate(rows) if r[0] in ("q1", "q2")][-2:]
    nrand = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    pick += [len(rows) - 1, len(rows) - 2] + random.sample(range(len(rows)), nrand)
    pick += sorted([i for i, r in enumerate(rows) if float(r[4]) > 10], key=lambda i: -float(rows[i][5]))[:2]
    pick += sorted([i for i, r in enumerate(rows) if float(r[4]) > 20], key=lambda i: -float(rows[i][5]))[:2]
    worst = 0
    for i in pick:
        k, a, c, x0, x1, v = rows[i]
        L0 = mp.mpf(a) if a != "-" else None; L1 = mp.mpf(c) if c != "-" else None
        m = my_cell(k, L0, L1, mp.mpf(x0), mp.mpf(x1))
        rel = (m - mp.mpf(v)) / mp.mpf(v); worst = max(worst, abs(rel))
        print("cell %4d %s L=[%s,%s] X=[%s,%s]: cert %s  indep %s  rel.diff %s" % (i, k, a[:7], c[:7], x0[:7], x1[:7], v, s(m, 8), s(rel, 3)), flush=True)
    print("max |rel diff| =", s(worst, 3))

def my_tail(XT, th=1):
    """re-derivation of DERIVATION 6.4 (own code)"""
    XT = mp.mpf(XT); dT = th * mp.exp(-b * XT)
    N0 = (1 + dT) ** B / B; K1 = (1 + dT) ** B * (1 + mp.exp(-B * XT / 2))
    assert 1 / N0 > S0 and 2 / (B * dT) > K1 ** 2 and rho0(K1 * B) < 0
    cth = mp.log((2 + dT) / (2 * pi * th)); lamX = b + max(cth, 0) / XT
    za_c = A57 * N0 + (2 * K1 + K1 * (mp.log(2 * pi) + mp.log(N0))) / pi - S0 / (pi * B)   # za <= za_c + K1 lam/pi
    zb_c = B57 / B + K1 * (CB / 2 * mp.log(1 + 2 * N0) + 1 / pi + CB / 2) + K1 * (mp.log(1 + 2 / UE) + dT / 2)
    T = mp.mpf(14)
    zz_c = K1 * (Nz(T) / 14 + CR1 / 14 + CR2 / (14 * mp.log(T)) + 1 / pi + CR1 / 28 + CR2 / (28 * mp.log(T)))
    eps = mp.exp((HALF - B) * XT)
    M = max(2 * (1 + dT) ** HALF * B, (1 + dT) ** HALF * (1 + mp.exp(-XT / 4)), ((1 + dT) ** mp.mpf(1.5) + 1) * (1 + mp.exp(-XT / 4)) / 2)
    fac = (1 + M * eps) / 2
    R = fac * ((za_c / XT + K1 * lamX / pi) / 4 + LOG9 * (za_c / XT + K1 * lamX / pi) / XT + 2 * N0 * mp.log(1 + XT) / XT ** 2
               + (max(zb_c, 0) + zz_c) / XT ** 2 + K1 * lamX ** 2 / pi + 2 * K1 * lamX / (pi * XT))
    xB = mp.exp(B * XT); lg = mp.log((1 + dT) / (1 - dT))
    R += (PSI_RS * eps + th * (1 + mp.exp(-XT / 2)) / 2) / XT ** 2
    R += (XT / 2 + lg) / ((1 - 4 / mp.exp(XT)) * xB * XT ** 2)
    R += ((XT / 2 + lg) / LOG2 + 1) / (4 * XT * xB)
    return R

def direct_point(L, X, th=1):
    """direct value of T1..T6 /(L+X)^2 at one (q = e^L, x = e^X), phi >= 1, delta = th e^{-bX}"""
    d = th * mp.exp(-b * X); al = mp.exp(-X / 2); xh = mp.exp((HALF - B) * X); xB = mp.exp(B * X)
    G, GL = trap(d, al)
    T1 = PSI_RS * xh; T2 = d * (mp.exp(b * X) + xh) / 2
    T3 = (zeta_sum(*G) + xh * zeta_sum(*GL)) / 2
    lg = mp.log((1 + d) / (1 - d))
    T5 = (X / 2 + lg) / (1 - (mp.exp(X / 2) * (1 - d)) ** -2) / xB
    T4 = 0; T6 = 0
    if L >= mp.log(3):
        Lq = max(L, LOG9)
        T4 = (Zof(zco(*G), G[0], Lq) + xh * Zof(zco(*GL), GL[0], Lq)) / 2
        T6 = ((L / LOG2) * (X / 2 + lg) + L) / xB
    return (T1 + T2 + T3 + T4 + T5 + T6) / (L + X) ** 2

def mode_ptail():
    for XT in (60, 100, 1000, 10**6):
        print("tail XT=%s: my re-derivation %s (cert: 0.01420, 0.01113, 0.007875, 0.007590)" % (XT, s(my_tail(XT), 6)), flush=True)
    lim = (b / (4 * pi) + b * b / pi) / 2
    print("limit (b/(4pi) + b^2/pi)/2 =", s(lim, 6))
    for X in (60, 61, 80, 150, 1000):
        X = mp.mpf(X); best = (0, None)
        for fr in (0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.87, 0.95, 1.0):
            L = X * fr
            if 0 < L < mp.log(3): L = mp.log(3)
            v = direct_point(L, X)
            if v > best[0]: best = (v, float(fr))
        print("direct X=%s: max over L/X in grid of the T1..T6 bound = %s at L/X=%s" % (X, s(best[0], 6), best[1]), flush=True)

def lp_val(L0, L1, kap, h):
    X = mp.log(kap * L0) / b
    if X < mp.log(2974): X = mp.log(2974)
    x = mp.exp(X)
    G = (1 / b, mp.mpf(1), ((1 + h) ** B + 1) / h); GL = (mp.mpf(2), mp.mpf(1), ((1 + h) ** HALF + 1) / h)
    Lq = max(L1, LOG9)
    Z = Zof(zco(*G), G[0], Lq); ZL = Zof(zco(*GL), GL[0], Lq)
    return X + mp.log(1 + h) - EUL + mp.mpf(0.2) / X + mp.mpf(0.2) / X ** 2 + (mp.exp(-b * X) * Z + mp.exp(-X / 2) * ZL) / 2 + 1 / (x - 1)

def mode_lprime():
    grid = [(mp.mpf(k) / 100, mp.mpf(hh) / 10) for k in (20, 25, 30, 40, 60, 100, 200) for hh in (3, 5, 7, 10)]
    L0 = mp.log(3); worst = None
    for i in range(6):
        L1 = L0 + max(mp.mpf(1) / 20, L0 / 100)
        v = min(lp_val(L0, L1, k, h) - INVB * mp.log(L0) for (k, h) in grid)
        print("Thm4 cell [%s,%s]: C = %s" % (s(L0, 5), s(L1, 5), s(v, 6)), flush=True)
        worst = v if worst is None or v > worst else worst
        L0 = L1
    print("Thm4 worst of first cells:", s(worst, 6), "(cert all cells 10.526)")

if __name__ == "__main__":
    t0 = time.time()
    md = sys.argv[1] if len(sys.argv) > 1 else "const"
    {"const": mode_const, "lemz": mode_lemz, "linnik": mode_linnik, "nonres": mode_nonres, "small": mode_small}.get(md, lambda: None)()
    if md in ("pnt", "ptail", "lprime"):
        globals()["mode_" + md]()
    print("[%s] %.1f s" % (md, time.time() - t0))
