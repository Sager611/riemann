#!/usr/bin/env python3
"""fals_zeta.py -- LiRangeAudit falsification for Li200k (zeta). One thread; Arb (python-flint) for bounds, numpy for zero sums.
T1 corner zeros: per-zero quartet inequality  2-2cosh(n l)cos(n phi) >= g(n theta) - q_n  at delta = a (edge), gamma just above T.
T2 lambda_n, n <= 200: Cauchy integral of d/dz log xi(1/(1-z)) on |z| = 0.85 (mpmath) vs Taylor series (Arb, as verify.py part3).
T3 small-T theorem (T = 100, 1000): lower bound L_n(T) = M_lo - Q_hi (Li200k formulas, HSW C3 = 9.4926, w-constant 0.1038+0.2573/log T)
   vs lambda_n (n <= 200); M_lo <= true M_n (Odlyzko zeros to 74920 + integral tail); Q_hi >= sum_{T<gamma<=Zmax} q_n;
   chain N_T; hypothetical configuration (zeros in (T, Zmax] moved to beta = 1/2 + a) lambda^hyp_n >= L_n(T) on a grid of n <= 1.3 N_T.
"""
import sys, time, math
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 200
PI = arb.pi(); EG = arb.const_euler()
def r(p, q=1): return arb(p) / q
C1, C2, C3 = r(1038, 10000), r(2573, 10000), r(94926, 10000) + r(7, 8)
def E(t): return C1 * t.log() + C2 * t.log().log() + C3
def th(t): return 2 * (1 / (2 * t)).atan()
def t_of(n): return (PI / (2 * arb(n))).cot() / 2
def K_consts(terms=80):
    K0 = arb(0); K1 = arb(0); lp = PI.log(); f = arb(1)
    for k in range(1, terms + 1):
        f = f * (2 * k - 1) * (2 * k); s = (-1) ** (k + 1); pw = PI ** (2 * k - 1)
        K0 += s * pw / (f * (2 * k - 1)); K1 += s * pw / f * (lp / (2 * k - 1) - arb(1) / (2 * k - 1) ** 2)
    tail = 2 * PI ** (2 * terms + 1) / (f * (2 * terms + 1) * (2 * terms + 2)) * (2 + lp)
    return K0 + arb(0, tail.upper()), K1 + arb(0, tail.upper())
K0, K1 = K_consts(); K1INF = PI / 2 * (1 - EG)
def main_term(n): n = arb(n); return n / PI * (K0 * ((n / (2 * PI)).log() - PI ** 2 / (4 * n * n)) - K1)
def M_lo(n1, n2):
    t1, t2 = t_of(n1), t_of(n2)
    corr = (C1 + C2 / t1.log()) * arb(n1) ** 2 / (2 * t1 * t1)
    a = main_term(n1) - 8 * E(t2) - corr
    if n1 >= 10 ** 4:
        tl = arb(n1).sqrt(); U1, U2 = arb(n1) * th(tl), arb(n2) * th(tl); n = arb(n1)
        b = n / PI * ((PI / 2 - 2 / U1) * ((n / (2 * PI)).log() - 1 / (4 * tl * tl)) - K1INF) - 8 * E(t2) * (U2 / PI + 3) - corr
        if b > a: return b
    return a
def q(n, a, t):
    x = arb(n) * a / (2 * t * t); return 4 * x.sinh() ** 2 + 4 * arb(n) * a * a / t ** 3
def Q_hi(n, a, T):
    cw = C1 + C2 / T.log()          # t E'(t) <= cw for t >= T  (verify.py uses 0.13, valid only for T >= 1.84e4)
    W = lambda t: t / (2 * PI) * (t / (2 * PI)).log() - t / (2 * PI) + cw * t.log()
    def tail(tK):
        xK = arb(n) * a / (2 * tK * tK); C = (xK.sinh() / xK) ** 2; L = (tK / (2 * PI)).log()
        i4 = (3 * L + 1) / (9 * tK ** 3) / (2 * PI) + cw / (4 * tK ** 4); i3 = (2 * L + 1) / (4 * tK ** 2) / (2 * PI) + cw / (3 * tK ** 3)
        return C * arb(n) ** 2 * a * a * i4 + 4 * arb(n) * a * a * i3
    X = arb(n) * a / (2 * T * T); bound = 2 * q(n, a, T) * E(T)
    if X < r(1, 1000): return bound + tail(T)
    Xf = float(X); xs = Xf; pts = [T]
    while xs > 1.0:
        xs = max(xs - 0.004, 1.0); pts.append(T * arb(Xf / xs).sqrt())
    t = pts[-1]
    while t < 100 * T: t = t * r(101, 100); pts.append(t)
    s = arb(0); Wp = W(pts[0])
    for k in range(len(pts) - 1):
        Wn = W(pts[k + 1]); s += q(n, a, pts[k]) * (Wn - Wp); Wp = Wn
    return bound + s + tail(pts[-1])
RATIOS = [2.0, 1.5, 1.2, 1.1, 1.05, 1.02, 1.01, 1.005, 1.002, 1.001, 1.0005, 1.0002, 1.0001, 1.00005, 1.00002, 1.00001]
def chain(a, T):
    n = 89
    while True:
        for rr in RATIOS:
            m = int(n * rr)
            if m > n and M_lo(n, m) - Q_hi(m, a, T) > 0: n = m; break
        else: return n

def g(x): return 2 - 2 * x.cos()
def T1():
    worst = None; cnt = 0
    for a in (r(1, 2), r(3, 8), r(37338196661320752, 10 ** 17), r(1, 8)):
        for T in (r(1415, 100), r(100), r(1000), r(10) ** 6, r(3000175332800)):
            for gam in (T + r(1, 10 ** 6), T * r(1001, 1000)):
                for d in (a, a * (1 - r(1, 10 ** 9)), a / 2):
                    if d >= r(1, 2): d = r(1, 2) - r(1, 10 ** 30)      # beta < 1 strictly
                    rho = acb(r(1, 2) + d, gam); z = 1 - 1 / rho; zs = 1 - 1 / (1 - rho.conjugate())
                    lz, lzs = z.log(), zs.log(); thg = th(gam)
                    nmax = 80 * float(gam) ** 2 / float(a)
                    ns = sorted(set(int(v) for v in np.geomspace(1, nmax, 160)) | set(int(nmax * f) + j for f in (0.5, 0.9, 1.0) for j in range(3)))
                    for n in ns:
                        lhs = (2 - (n * lz).exp() - (n * lzs).exp()).real        # per-zero share of the quartet
                        qq = q(n, a, gam); rhs = g(n * thg) - qq
                        dlt = lhs - rhs; cnt += 1
                        if not (dlt > 0):
                            print("T1 VIOLATION/UNDECIDED a=%s T=%s gam=%s d=%s n=%d lhs-rhs=%s" % (a.str(5), T.str(8), gam.str(15), d.str(5), n, dlt.str(5))); return
                        rel = float(dlt / qq) if qq > 0 else 1.0
                        if worst is None or rel < worst[0]: worst = (rel, float(a), float(T), float(gam), float(d), n)
    print("T1 corner zeros: %d certified cases, all lhs > rhs; min slack/q_n = %.3g at a=%.4g T=%.4g gamma=%.8g delta=%.4g n=%d" % ((cnt,) + worst))
def taylor_lambda(NMAX):
    from flint import arb_series
    old = ctx.prec; ctx.prec = int(1.1 * NMAX) + 120; ctx.cap = NMAX + 1
    x = arb_series([0, 1]); one = arb_series([1])
    f = (-arb_series.zeta(x)).log() + (one + x / 2).lgamma() + (one - x).log() - x * arb.pi().log() / 2
    co = f.coeffs(); lam = [None]
    for n in range(1, NMAX + 1):
        s = arb(0); c = arb(1)
        for k in range(1, n + 1):
            if k > 1: c = c * (n - k + 1) / (k - 1)
            s += (-1) ** k * c * co[k]
        lam.append(n * s)
    ctx.prec = old
    return lam
def cauchy_lambda(NMAX, M=512, rad=0.85):
    import mpmath as mp
    mp.mp.dps = 45
    def f(z):
        s = 1 / (1 - z)
        return s * s * (1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, 1, 1) / mp.zeta(s))
    vals = {}
    for k in range(M // 2 + 1):
        z = rad * mp.expjpi(2 * mp.mpf(k) / M); vals[k] = (z, f(z))
    lam = [None]
    for n in range(1, NMAX + 1):
        s = 0
        for k in range(M):
            z, fv = vals[k] if k <= M // 2 else (mp.conj(vals[M - k][0]), mp.conj(vals[M - k][1]))
            s += fv * z ** (-(n - 1))
        lam.append(mp.re(s) / M)
    return lam


def cauchy_lambda(NMAX, M=512, rad=0.85):         # (overrides the draft above) incremental powers
    import mpmath as mp
    mp.mp.dps = 45
    def f(z):
        s = 1 / (1 - z)
        return s * s * (1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, 1, 1) / mp.zeta(s))
    zk = []; fk = []
    for k in range(M // 2 + 1):
        z = rad * mp.expjpi(2 * mp.mpf(k) / M); zk.append(z); fk.append(f(z))
    for k in range(M // 2 + 1, M):
        zk.append(mp.conj(zk[M - k])); fk.append(mp.conj(fk[M - k]))
    zinv = [1 / z for z in zk]; pw = [mp.mpf(1)] * M; lam = [None]
    for n in range(1, NMAX + 1):
        lam.append(mp.re(mp.fsum(fk[k] * pw[k] for k in range(M))) / M)
        pw = [pw[k] * zinv[k] for k in range(M)]
    return lam
def main():
    t_start = time.time()
    T1(); print("  (%.0f s)" % (time.time() - t_start), flush=True)
    lt = taylor_lambda(200); lc = cauchy_lambda(200)
    dmax = max(abs(float(lt[n]) - float(lc[n])) / abs(float(lt[n])) for n in range(1, 201))
    print("T2 lambda_n n<=200: Cauchy (|z|=0.85, 512 pts, 45 dig) vs Taylor (Arb): max rel diff %.2e; lambda_1 = %s, lambda_100 = %.10f, lambda_200 = %.10f; min_n lambda_n = %.6g"
          % (dmax, lt[1].str(15), float(lt[100]), float(lt[200]), min(float(lt[n]) for n in range(1, 201))), flush=True)
    import mpmath as mp
    Z = np.loadtxt(sys.argv[1]); Zmax = float(Z[-1]); thZ = 2 * np.arctan(1 / (2 * Z))
    EZ = float(E(arb(Zmax)))
    cache = {}
    def tail_int(n):
        if n in cache: return cache[n]
        mp.mp.dps = 20; U = n * 2 * math.atan(1 / (2 * Zmax))
        h = lambda u: (1 - mp.cos(u)) / mp.pi * mp.log(mp.cot(u / (2 * n)) / 2 / (2 * mp.pi)) / (4 * n) / mp.sin(u / (2 * n)) ** 2
        pts = [0] + [k * math.pi for k in range(1, int(U / math.pi) + 1)] + [U]
        v = float(mp.quad(h, pts)); err = 8 * EZ * (U / math.pi + 2) + 0.13 * n * n / (2 * Zmax ** 2)
        cache[n] = (v, err); return cache[n]
    def M_true(n):
        v, err = tail_int(n); return float(np.sum(2 - 2 * np.cos(n * thZ))) + v, err
    def Q_true(n, a, T):
        y = Z[Z > T]; x = n * a / (2 * y * y); return float(np.sum(4 * np.sinh(x) ** 2 + 4 * n * a * a / y ** 3))
    def lam_hyp(n, a, T):
        lo = Z[Z <= T]; hi = Z[Z > T]; m = len(hi) // 2 * 2; y = hi[0:m:2]
        d = min(a, 0.5 - 1e-15)
        l = 0.5 * np.log(((0.5 + d) ** 2 + y * y) / ((0.5 - d) ** 2 + y * y)); ph = np.arctan((0.5 - d) / y) + np.arctan((0.5 + d) / y)
        s = float(np.sum(2 - 2 * np.cos(n * 2 * np.arctan(1 / (2 * lo))))) + float(np.sum(4 - 4 * np.cosh(n * l) * np.cos(n * ph)))
        s += float(np.sum(2 - 2 * np.cos(n * 2 * np.arctan(1 / (2 * hi[m:])))))
        v, err = tail_int(n); return s + v, err
    for Tv in (100, 1000):
        T = arb(Tv)
        for a, an in ((r(1, 2), 0.5), (r(3, 8), 0.375)):
            t1 = time.time(); NT = chain(a, T)
            bad = []
            for n in range(89, 201):
                if not (float((M_lo(n, n) - Q_hi(n, a, T)).upper()) <= float(lt[n])): bad.append(n)
            grid = sorted(set(int(v) for v in np.geomspace(89, 3 * NT, 110)) | {NT, NT + 1})
            fM = fQ = fL = 0; minhyp = None; first_neg = None
            for n in grid:
                Lb = M_lo(n, n) - Q_hi(n, a, T); Lf = float(Lb.upper())
                Mt, err = M_true(n)
                if float(M_lo(n, n).lower()) > Mt + err: fM += 1
                if float(Q_hi(n, a, T).upper()) < Q_true(n, an, T): fQ += 1
                lh, err2 = lam_hyp(n, an, T)
                if lh + err2 < Lf: fL += 1
                if n <= NT and (minhyp is None or (lh - err2) / Mt < minhyp[0]): minhyp = ((lh - err2) / Mt, n)
                if lh + err2 < 0 and first_neg is None: first_neg = n
            print("T3 T=%d a=%s: chain N_T = %d = %.3f T^2; L_n <= lambda_n fails for n in 89..200: %s; on %d grid n in [89, 3N_T]: M_lo>M_true %d, Q_hi<Q_true %d, lambda^hyp<L_n %d; min_{n<=N_T} (lambda^hyp-err)/M_n = %.3g at n=%d; first grid n with lambda^hyp+err<0: %s (%.0f s)"
                  % (Tv, a.str(3), NT, NT / Tv ** 2, bad or 'none', len(grid), fM, fQ, fL, minhyp[0], minhyp[1], first_neg, time.time() - t1), flush=True)
    print("total %.0f s" % (time.time() - t_start))
if __name__ == "__main__":
    main()

