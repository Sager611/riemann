# t3b_universal.py — explore6 T3: methods that use only "g multiplicative, |g| <= 1" (Katai/BSZ, Daboussi,
# Montgomery-Vaughan) cannot give power savings even against a badly approximable rotation.
# g(p) = e(-p*phi) on primes, g(p^k) = 0 (k >= 2), so g(n) = mu^2(n) e(-phi * sopf(n)); then
# T(x) = sum_{n<=x} g(n) e(n*phi) contains pi(x) exactly (phase 0 at primes) and |T(x)| ~ x/log x.
# Compare: sum_{n<=x} mu(n) e(n*phi) ~ sqrt(x) (Test A), proven x^{4/5+eps} (Vaughan-type minor-arc bound).
# Memory: int8/int32 sieve arrays, chunked sums (peak well under 1 GB).
import numpy as np, math, time
t0 = time.time()
X = 10**7
mu = np.ones(X + 1, np.int8); prod = np.ones(X + 1, np.int32); sopf = np.zeros(X + 1, np.int32)
R = int(math.isqrt(X)) + 1
sv = np.ones(R + 1, bool); sv[:2] = False
for i in range(2, int(R**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
for p in np.nonzero(sv)[0].tolist():
    mu[::p] *= -1; prod[::p] *= p; sopf[::p] += p; mu[::p*p] = 0
phi = (math.sqrt(5) - 1) / 2
checks = {10**4, 10**5, 10**6, 10**7}
T = 0j; Smu = 0j; pis = 0; res = {}; Msum = 0
CH = 10**6
for lo in range(0, X + 1, CH):
    hi = min(lo + CH, X + 1)
    n = np.arange(lo, hi, dtype=np.int64)
    pr = prod[lo:hi].astype(np.int64)
    q = n // np.maximum(pr, 1)
    big = q > 1
    m = mu[lo:hi].copy(); m[big] *= -1
    so = sopf[lo:hi].astype(np.int64); so[big] += q[big]
    if lo == 0: m[0] = 0
    isp = (so == n) & (n >= 2)
    sq = m != 0
    g = np.where(sq, np.exp(2j * math.pi * np.mod((n - so).astype(np.float64) * phi, 1.0)), 0)
    f = m.astype(np.float64) * np.exp(2j * math.pi * np.mod(n.astype(np.float64) * phi, 1.0))
    cT = np.cumsum(g) + T; cS = np.cumsum(f) + Smu; cp = np.cumsum(isp) + pis
    for x in checks:
        if lo <= x < hi:
            res[x] = (int(cp[x - lo]), cT[x - lo], cS[x - lo])
    T, Smu, pis = cT[-1], cS[-1], int(cp[-1])
    if lo == 0:
        Msum = int(m.astype(np.int64).sum())          # sum_{n < 10^6} mu(n)
    if lo == CH:
        assert Msum + int(m[0]) == 212                  # M(10^6) = 212
print("x         pi(x)     |T(x)| (g universal)   |T|/pi(x)   |T - pi|/sqrt(x)   |sum mu e(n phi)|/sqrt(x)")
for x in sorted(res):
    p_, t_, s_ = res[x]
    print(f"{x:<9d} {p_:<9d} {abs(t_):<22.1f} {abs(t_)/p_:.4f}      {abs(t_-p_)/math.sqrt(x):.3f}             {abs(s_)/math.sqrt(x):.3f}")
print(f"seconds {time.time()-t0:.1f}")
