#!/usr/bin/env python3
"""Numerical check of Lemma prelim:prime-estimates (TeX l.308-407) for x <= 1e7:
psi(x) <= log4 * x, M(x) <= log x + log 4, A(n) two-sided bounds (prelim:von-mangoldt-bounds),
and the interval bound M(x)-M(y) > log(x/y) - 2 for 100 <= y <= x (worst case reported).
Memory-lean: only prime / prime-power arrays (step functions are checked at their jump points)."""
import numpy as np, math
N = 10**7
s = np.ones(N + 1, dtype=np.bool_); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.nonzero(s)[0].astype(np.int64); del s
lp = np.log(P.astype(np.float64))
# prime powers p^k <= N with Lambda = log p, sorted
pp = [P]; lw = [lp]
for k in range(2, 24):
    m = P[P.astype(np.float64) ** k <= N]
    if len(m) == 0: break
    pp.append(m ** k); lw.append(np.log(m.astype(np.float64)))
pp = np.concatenate(pp); lw = np.concatenate(lw); o = np.argsort(pp); pp = pp[o]; lw = lw[o]
psi = np.cumsum(lw)                       # psi at each prime power (step fn, max of psi(n)/n at jumps)
r = psi / pp
i = int(np.argmax(r[pp >= 2]))
print("max psi(n)/n, 2<=n<=1e7:", float(r.max()), "at n =", int(pp[i]), " vs log4 =", math.log(4))
M = np.cumsum(lp / P)
print("max M(p)-log p (sup of M(x)-log x):", float((M - lp).max()), " vs log4 =", math.log(4))
A = np.cumsum(lw / pp)                    # A at each prime power
ppf = pp.astype(np.float64)
hi = np.log(ppf) - 1 + math.log(4) + (1 + np.log(ppf)) / ppf   # hi increasing: min of hi-A at jumps
print("min over n of upper-A(n):", float((hi - A).min()))
nr = np.append(ppf[1:] - 1, N)            # right ends of constancy intervals; lo increasing
lo = np.log(nr) - 1 + 1 / nr
print("min over n of A(n)-lower:", float((A - lo).min()))
# interval bound: E(t)=M(t)-log t; sup on [p_k,p_{k+1}) at p_k, inf -> p_{k+1}
E_sup = M - lp
E_inf = M[:-1] - lp[1:]
i100 = int(np.searchsorted(P, 100))
E100 = M[i100 - 1] - math.log(100)
runmax = np.maximum.accumulate(np.concatenate(([E100], E_sup[i100:])))
cand0 = (M[i100 - 1] - math.log(101)) - E100
cand = E_inf[i100:] - runmax[1:len(E_inf[i100:]) + 1]
worst = min(cand0, float(cand.min()))
k = int(np.argmin(cand))
print("min over 100<=y<=x<=1e7 of M(x)-M(y)-log(x/y):", worst,
      " (claim > -2); attained near x =", int(P[i100 + k + 1]))
print("sum_{odd p<=1e7} log p/p^2 =", float(np.sum(lp[1:] / P[1:].astype(np.float64)**2)), "(claim < 7/20)")
print("V = sum log p/(p(p-1)) =", float(np.sum(lp / (P * (P - 1.0)))), "(claim < 33/40)")

