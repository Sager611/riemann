#!/usr/bin/env python3
"""Data tests of Prop. real:hyperbola and Lemma real:large-prime (TeX l.846-1134), arXiv 2610.30476.
No real character has a known real zero, so we test the proof's identities WITHOUT assuming L(beta)=0:
running the paper's proof verbatim (it uses L(beta)=0 only at (real:first-half) and (real:paired-tails)) gives
  H_beta(t) = L(1)t^eps/eps + L(beta)*G(t) + E(t),  G(t) = sum_{d<=sqrt t} d^-beta - t^{eps/2}/eps,
with the SAME claimed bound |E(t)| <= (q+1+log(t)/2) t^{1/2-beta}; likewise each of the three pieces.
For the Lemma: N(xy)-N(y) = L(1)kappa + L(beta)*Gterm + R with |R| <= B(xy)log(xy)+B(y)log y+2B(y)+2y^-1/2."""
import numpy as np, math, mpmath as mp, sys
mp.mp.dps = 40
T = int(sys.argv[1]) if len(sys.argv) > 1 else 2 * 10**6
def kron(D, n):
    if n == 0: return 0
    res = 1; m = n
    while m % 2 == 0:
        m //= 2
        if D % 2 == 0: return 0
        res *= 1 if D % 8 in (1, 7) else -1
    p = 3
    while p * p <= m:
        while m % p == 0:
            m //= p; r = pow(D % p, (p - 1) // 2, p)
            if D % p == 0: return 0
            res *= 1 if r == 1 else -1
        p += 2
    if m > 1:
        if D % m == 0: return 0
        res *= 1 if pow(D % m, (m - 1) // 2, m) == 1 else -1
    return res
def Lval(chi_per, q, s):
    if s == 1:
        return -mp.fsum(chi_per[a] * mp.digamma(mp.mpf(a) / q) for a in range(1, q)) / q
    return mp.fsum(chi_per[a] * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q)) / mp.mpf(q) ** s
worst_all = {}
for D in [-3, -4, 5, -7, 8, -8, 12, 13, -163]:
    q = abs(D)
    chi_per = [kron(D, a) for a in range(q)]
    n = np.arange(T + 1, dtype=np.float64)
    chi = np.array([chi_per[i % q] for i in range(T + 1)], dtype=np.float64); chi[0] = 0
    a = np.zeros(T + 1)
    for dd in range(1, T + 1):
        c = chi[dd]
        if c: a[dd::dd] += c
    a[0] = 0
    # a = 1*chi: sum over divisors d of chi(d): a[n] += chi(d) for d | n  (loop above)
    assert a[1] == 1 and a.min() >= 0
    L1 = Lval(chi_per, q, 1)
    for beta in [0.5, 0.6, 0.75, 0.9, 0.99, 0.999]:
        eps = 1 - beta
        Lb = Lval(chi_per, q, mp.mpf(beta))
        L1f, Lbf = float(L1), float(Lb)
        w = np.zeros(T + 1); w[1:] = n[1:] ** (-beta)
        Hb = np.cumsum(a * w)
        K = math.isqrt(T) + 2
        Z = np.cumsum(np.concatenate(([0.0], np.arange(1, K + 1, dtype=np.float64) ** (-beta))))
        nn = n[1:]; isq = np.array([math.isqrt(int(v)) for v in nn]) if T <= 10**5 else np.floor(np.sqrt(nn) + 1e-9).astype(np.int64)
        isq = np.where((isq + 1) ** 2 <= nn, isq + 1, isq); isq = np.where(isq ** 2 > nn, isq - 1, isq)
        Gt = Z[isq] - nn ** (eps / 2) / eps
        E = Hb[1:] - L1f * nn ** eps / eps - Lbf * Gt
        bound = (q + 1 + 0.5 * np.log(nn)) * nn ** (0.5 - beta)
        r_right = np.abs(E) / bound
        # left limits t -> (n+1)^-: H_beta(n), floor(sqrt t) = isqrt(n)
        t2 = nn + 1
        Gl = Z[isq] - t2 ** (eps / 2) / eps
        El = Hb[1:] - L1f * t2 ** eps / eps - Lbf * Gl
        r_left = np.abs(El) / ((q + 1 + 0.5 * np.log(t2)) * t2 ** (0.5 - beta))
        rmax = max(float(r_right.max()), float(r_left.max()))
        # three pieces at sampled t (O(sqrt t) each)
        S1 = np.cumsum(chi / np.maximum(n, 1)); Sb = np.cumsum(chi * w)
        Zfull = np.cumsum(w)
        pr = [0.0, 0.0, 0.0, 0.0]
        for t in np.unique(np.concatenate((np.arange(1, 400), np.geomspace(400, T // 2, 300).astype(np.int64)))):
            t = float(t); u = math.sqrt(t); U = int(math.floor(u + 1e-12))
            if (U + 1) ** 2 <= t: U += 1
            ds = np.arange(1, U + 1)
            Rb_td = Lbf - Sb[(t // ds).astype(np.int64)]
            F1 = -float(np.sum(ds ** (-beta) * Rb_td))
            ms = ds
            second = float(np.sum(chi[ms] * ms ** (-beta) * (Zfull[(t // ms).astype(np.int64)] - Zfull[U])))
            R1u = L1f - S1[U]; Rbu = Lbf - Sb[U]
            main2 = (t ** eps / eps) * (L1f - R1u) - (u ** eps / eps) * (Lbf - Rbu)
            E0 = second - main2
            E1 = -(t ** eps / eps) * R1u + (u ** eps / eps) * Rbu
            tb = t ** (0.5 - beta)
            pr[0] = max(pr[0], abs(F1) / (q / 2 * tb))
            pr[1] = max(pr[1], abs(E0) / ((1 + 0.5 * math.log(t)) * tb))
            pr[2] = max(pr[2], abs(E1) / (q / 2 * tb))
            # identity check: H_beta(t) = L(beta) Z(u) + F1 + second
            ident = Hb[int(t)] - (Lbf * Zfull[U] + F1 + second)
            pr[3] = max(pr[3], abs(ident) / max(1.0, abs(Hb[int(t)])))
        worst_all[(D, beta)] = (rmax, pr)
        print(f"D={D:5d} beta={beta:<5}  max|E|/bound={rmax:.4f}  pieces/bounds: first-half={pr[0]:.3f} E0={pr[1]:.3f} E1={pr[2]:.3f}  identity err={pr[3]:.1e}")
    # ---- Lemma real:large-prime, analytic part (generalised), at beta=0.99 and 0.999 where eps log(xy)<=1
    logn = np.zeros(T + 1); logn[1:] = np.log(n[1:])
    Ncum = np.cumsum(a * logn / np.maximum(n, 1))
    Hone = np.cumsum(a / np.maximum(n, 1))
    for beta in [0.99, 0.999]:
        eps = 1 - beta; Lbf = float(Lval(chi_per, q, mp.mpf(beta))); L1f = float(L1)
        Zd = np.cumsum(np.concatenate(([0.0], np.arange(1, math.isqrt(T) + 3, dtype=np.float64) ** (-beta))))
        phi = lambda t: math.log(t) * t ** (-eps)
        def Gint(y, X):   # int_y^X G(t) phi'(t) dt, G = Z(floor sqrt t) - t^{eps/2}/eps
            tot = 0.0; k = math.isqrt(int(y))
            while k * k <= X:
                lo = max(k * k, y); hi = min((k + 1) ** 2, X)
                if hi > lo: tot += Zd[k] * (phi(hi) - phi(lo))
                k += 1
            Fm = lambda t: 2 * t ** (-eps / 2) * (1 / eps ** 2 + math.log(t) / eps)   # antiderivative of t^{eps/2}phi'/eps
            return tot - (Fm(X) - Fm(y))
        Gf = lambda t: Zd[math.isqrt(int(t))] - t ** (eps / 2) / eps
        Bf = lambda t: (q + 1 + 0.5 * math.log(t)) * t ** -0.5
        worstR = 0.0; worstJ = -1e9; cnt = 0
        for y in [100, 300, 1000, 3000]:
            for x in [y, 2 * y, 10 * y, T // y]:
                if x < y or x * y > T or eps * math.log(x * y) > 1: continue
                X = x * y
                kap = math.log(x) * math.log(y) + 0.5 * math.log(x) ** 2
                Gterm = phi(X) * Gf(X) - phi(y) * Gf(y) - Gint(y, X)
                R = (Ncum[X] - Ncum[y]) - L1f * kap - Lbf * Gterm
                Rb = Bf(X) * math.log(X) + Bf(y) * math.log(y) + 2 * Bf(y) + 2 * y ** -0.5
                worstR = max(worstR, abs(R) / Rb)
                # weighted injection (unconditional) and the generalised final inequality
                ps = [p for p in range(y + 1, x + 1) if a[p] > 0 and all(p % k for k in range(2, math.isqrt(p) + 1))]
                J = sum(a[p] * math.log(p) / p for p in ps)
                inj = Hone[y] * J <= Ncum[X] - Ncum[y] + 1e-9
                Ecal = Bf(y) * (eps * kap + math.log(y)) + Bf(X) * math.log(X) + 2 * Bf(y) + 2 * y ** -0.5
                extra = Lbf * (Gterm - eps * kap * y ** (-eps) * Gf(y))
                rhs = eps * kap + (Ecal + extra) / Hone[y]
                worstJ = max(worstJ, J - rhs); cnt += 1
                assert inj, (D, beta, x, y)
        print(f"   Lemma (generalised), beta={beta}: {cnt} (x,y) pairs: max|R|/bound={worstR:.4f}; max(J - rhs)={worstJ:.3e} (must be <=0); weighted injection held")
print("OVERALL max |E|/bound =", max(v[0] for v in worst_all.values()),
      " max piece ratio =", max(max(v[1][:3]) for v in worst_all.values()))

