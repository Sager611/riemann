#!/usr/bin/env python3
"""tail_margin.py (ppraudit) -- independent (mpmath, own kernel/Lemma-Z/sup code from indep_kappa.py) recomputation of
(a) the tail final inequality (5.6) at c_T = 391/200 with the certificate's tau, K_A, K_B, c_sigma;
(b) float re-integration of I, J, K3 (own psi, midpoint rule on 2e5 cells) and the u >= 60 bounds A_rel, B_rel;
(c) Corollary E binding cases: w = 27 at L = 94.4782 (g* < p), w = 52 at L = 224.609 (g* < sqrt p - 2), own least kappa."""
import sys, math
sys.argv = [sys.argv[0]]
import indep_kappa as IK
M = IK.M; M.mp.prec = 256
# (a)
tau = M.mpf("0.101433") - M.mpf("1e-6") * 0 ; tau = M.mpf("0.1014329")   # certificate lower endpoint (0.101433 +- 9.69e-8)
KA, KB, csig = M.mpf("0.0160555"), M.mpf("0.00906309"), M.mpf("3.2795975")
P = IK.PR[:10001]; L1 = M.fsum(M.log(q) for q in P); u1 = M.log(P[-1])       # log(P+1) - theta < 1e-40000
s, d = M.mpf(23) / 10, M.mpf(1) / 10; beta = 1 + s + d
def tailmargin(cT):
    X = beta * M.log(cT * L1) * IK.IB; Y = X - IK.ell; eY = M.e * Y
    k1, k2 = IK.n1 + IK.n0 / eY, IK.n2 + 2 * IK.n1 / eY + (2 / eY) ** 2 * IK.n0
    k1L, k2L = IK.m1 + IK.m0 / eY, IK.m2 + 2 * IK.m1 / eY + (2 / eY) ** 2 * IK.m0
    ZoL = IK.supZL(IK.zco(IK.n0, k1, k2), IK.n0, L1); ZLoL = IK.supZL(IK.zco(IK.m0, k1L, k2L), IK.m0, L1)
    xB, xh = M.e ** (IK.B * X), M.e ** ((M.mpf(1) / 2 - IK.B) * X)
    z1, z1L = IK.SIG * k2 / 2, IK.SIG * k2L / 2
    triv = IK.W0 / ((1 - M.e ** (-2 * Y)) * xB * L1)
    sub = (z1 + xh * z1L) / L1 + triv + IK.winf * (7 / L1 + 1) / xB
    Ro = ZoL / 2 + xh * ZLoL / 2 + triv
    kz = 1 - 1 / u1; ell1 = M.log(L1 / kz); mdel = L1 ** d / ell1
    taup = tau * M.e ** (-M.euler) * (1 - 1 / (2 * u1 * u1))
    lhs = taup * (IK.W1 * cT ** beta * mdel - sub / (L1 ** s * ell1))
    rhs = Ro * (M.mpf(2) ** 15 / L1 ** s + csig ** 2 * (KA + KB) / kz ** s)
    return lhs - rhs, lhs, rhs
mg, lh, rh = tailmargin(M.mpf(391) / 200)
print("(a) tail (5.6) at c_T = 1.955: margin %s (lhs %s, rhs %s; certificate margin 0.0004080)" % (M.nstr(mg, 8), M.nstr(lh, 8), M.nstr(rh, 8)))
mg2, _, _ = tailmargin(M.mpf(1954) / 1000)
print("    at c_T = 1.954: margin %s" % M.nstr(mg2, 6))
# (b) own psi (float), integrals
EG = math.exp(0.5772156649015329); LQ = math.log(1e5); epsM = 1 / (2 * LQ ** 2); etaB = 1 / LQ ** 2 + 1e-5
small = [p for p in IK.PR if p <= 283]
def ein(c): return sum(c ** k / (k * math.factorial(k)) for k in range(1, 60))
RB = {c: sum((p ** (c / LQ) - 1) / p for p in small) + ein(c) + (math.exp(c) - 1) * (1 / (2 * LQ ** 2) + 1 / (2 * math.log(288) ** 2)) for c in (0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6)}
def psi(a):
    v = [EG / ((1 - epsM) * a) - 1] if a <= 1 else [EG / (1 - epsM) - 1]
    if a > 1:
        am = min(a, 2.609); h = 2 * am - am * math.log(am) - 1 - etaB * (am - 1) - (math.log(am) + etaB) / LQ
        if h > 0: v.append(EG / ((1 - epsM) * h) - 1)
    for c, R in RB.items():
        r = math.exp(-c * a + R)
        if r < 1: v.append(r / (1 - r))
    return min(v)
Psi = lambda t: psi((2.3 / t - 1) / 2)
n = 200000; I = J = 0.0; K3 = 0.0
for i in range(n):
    t = (i + 0.5) / n; ps = Psi(t); I += ps / t ** 2 / n; J += ps / t ** 4 / n; K3 = max(K3, ps / t ** 3)
print("(b) own midpoint integrals: I = %.6f (cert upper 0.882775), J = %.5f (cert 3.27487), K3 = %.5f (cert 4.72882), Psi(1) = %.5f (cert 1.75049)" % (I, J, K3, Psi(1.0)))
epsV = (1 + 1 / (2 * LQ ** 2)) / (1 - 1e-5)
Brel = epsV * (I + (2 * Psi(1.0) + 4 * J + K3) / 7200) / (1 - 1 / 7200)
print("    B_rel(u = 60) with own integrals = %.6f (cert 0.889211)" % Brel)
# (c) Corollary E binding cases
for w, Lc, fac, sl in ((27, M.mpf("94.4782"), 1, 0), (52, M.mpf("224.609"), M.mpf(1) / 2, M.mpf(1) / 1000)):
    y = IK.tab[w][0]; Tn, C, _ = IK.rosser(w, y); T = M.mpf(Tn) / IK.ONE
    lo_, hi_ = M.mpf(1), M.mpf(10) ** 9
    for _ in range(60):
        mid = M.sqrt(lo_ * hi_); mm = IK.margin(mid, Lc, T, C)
        if mm is not None and mm > 0: hi_ = mid
        else: lo_ = mid
    lhs = M.log(hi_ * Lc) * IK.IB + sl
    print("(c) w = %d at L = %s: own least kappa %s; (1/b) log(kappa L) + slack = %s vs %s*L = %s -> %s" % (w, Lc, M.nstr(hi_, 6), M.nstr(lhs, 8), fac, M.nstr(fac * Lc, 8), "holds" if lhs < fac * Lc else "FAILS"))
