# probes (c), (d): exponent bookkeeping for ThinH and the scale/uniformity of Y1 in PairEx (pure python, log-space)
import math
print('(c) ThinH: classical saving c*log y/(lam + sqrt(log y)) at y = e^{x^{eps/3}}/3, lam = 2 x^{eps/5} log^2 x, vs required x^{eps/8} (archive c = 1)')
def need_L(eps, c, C=10.0):
    # smallest L = log x with  c*ly/(lam+sqrt(ly)) >= x^{eps/8} + 2 log(ly+lam) + log C, in log space
    # ly ~ x^{eps/3} (the -log 3 is negligible), lam = 2 x^{eps/5} L^2
    def lae(a, b): return max(a, b) + math.log1p(math.exp(-abs(a - b)))
    def f(L):
        lly = eps / 3 * L; llam = math.log(2) + eps / 5 * L + 2 * math.log(L)
        A = math.log(c) + lly - lae(llam, lly / 2)
        extra = 2 * math.log(math.exp(min(lly, 700)) + math.exp(min(llam, 700))) + math.log(C) if max(lly, llam) < 700 else 2 * max(lly, llam) + math.log(C)
        B = lae(eps / 8 * L, math.log(max(extra, 1e-300)))
        return A - B
    lo, hi = 20.0, 20.0
    while f(hi) < 0: hi *= 1.5
    lo = hi / 1.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0: lo = mid
        else: hi = mid
    return hi
for eps in (1/52, 1/104, 1/520):
    for c in (0.1, 0.02):
        L = need_L(eps, c)
        print('  eps=1/%d c_dlVP=%.2f: holds from log x ~ %.3g (x ~ 10^%.3g); exponent margin 2eps/15 - eps/8 = eps/120 = %.2e' % (round(1/eps), c, L, L / math.log(10), eps / 120))
print('  Siegel zero instead of RealBox: 1-beta >= c(d) N^{-d}: (1-beta) log y <= c * exp(-2d x^{eps/5} log^2 x) * x^{eps/3} -> 0: no saving; RealBox is necessary for chars with q|N or d|N.')
print()
print('(c2) j=1 box height vs classical zero density (Montgomery 1971, (Q^2 T)^{2(1-s)/s} (log QT)^14, s = 1-1/400) and the archive count C Q^{1/4}, family |P_d| >= e^{0.9 z\'}')
for Ylog10 in (600, 1000, 3000):
    lY = Ylog10 * math.log(10); eps = 1/52
    z = math.exp(eps / 5 * lY); zp = z * lY ** 2; M = 1
    lQ = math.log(8 * M) + z + zp        # log height(M, eps, Y, 1) = log of the conductor bound
    lmont = 3 * lQ * (2 * (1/400) / (1 - 1/400)) + 14 * math.log(2 * lQ)
    print('  Y=10^%d: z=%.4g, log d ~ z log^2 Y=%.4g; log height_j=1 = %.4g; log Montgomery count = %.4g (=%.4f log Q); log(C Q^{1/4}) = %.4g + log C; log|P_d| >= %.4g'
          % (Ylog10, z, zp, lQ, lmont, lmont / lQ, lQ / 4, 0.9 * zp))
print()
print('(d) Scale of Y1 (PairEx) from the Markov step alone: (32 C/tau) Y^{-tau} <= eta/(2R), tau = eps/2, eta = 1/2, R = #chars mod 8M (x2 for a in {0,1})')
for eps in (1/52, 1/104):
    tau = eps / 2
    for M in (1, 4, 100):
        R = 2 * sum(1 for a in range(1, 8 * M + 1) if math.gcd(a, 8 * M) == 1)
        for C in (1.5, 2.5):
            lY = math.log(32 * C / tau * 2 * R / 0.5) / tau
            print('  eps=1/%d M=%d R=%d C=%.1f: Y1 >= e^%.4g = 10^%.4g' % (round(1/eps), M, R, C, lY, lY / math.log(10)))
print('  pointwise margin: y^{tau} >= r sqrt(log y) (r = 3, typical max ratio from probe (a)): log y >= %.4g' % (lambda tau: (lambda f: f(f, 300.0, 60))(lambda g, L, k: L if k == 0 else g(g, (math.log(3) + 0.5 * math.log(L)) / tau, k - 1)))(1/104))

print()
print('(e) principal case of gChi: G(y) = sum_{n<=y} Lambda(n)(y/n - 1) - y(log y - 1 - gamma) = -y int_y^inf (psi(u)-u)u^{-2} du; check G(y)/y -> 0 and |G(y)| vs sqrt(y) (numpy sieve to 1e7)')
import numpy as np
Zs = 10**7
sv = np.ones(Zs + 1, dtype=bool); sv[:2] = False
for i in range(2, int(Zs ** 0.5) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0]
ns, ls = [P.astype(np.float64)], [np.log(P.astype(np.float64))]
pk = P.astype(np.int64) ** 2
lp = np.log(P.astype(np.float64)); bp = P.astype(np.int64)
while True:
    m = pk <= Zs
    if not m.any(): break
    ns.append(pk[m].astype(np.float64)); ls.append(lp[m])
    pk = pk[m] * bp[m]; lp = lp[m]; bp = bp[m]
N = np.concatenate(ns); L = np.concatenate(ls); o = np.argsort(N); N = N[o]; L = L[o]
cL = np.cumsum(L); cLn = np.cumsum(L / N)
g = 0.5772156649015329
for y in (1e3, 1e4, 1e5, 1e6, 1e7):
    k = np.searchsorted(N, y, side='right')
    G = y * cLn[k-1] - cL[k-1] - y * (math.log(y) - 1 - g)
    print('  y=%.0e: G(y)=%.3f  G/y=%.2e  G/sqrt(y)=%.3f' % (y, G, G / y, G / math.sqrt(y)))

