#!/usr/bin/env python3
# b3check.py -- r2aaudit: numerical sanity check of Ramare's b3 = 0.605 inside KLN (their (3.7), [Ram, Lemma 5.6]):
#   sum_{X < n < 5X} lambda_X(n)^2 / n^2 <= b3 / X,  lambda_X(n) = sum_{d | n, d <= X} mu(d)   (stated for X >= 1e9).
# For X < n < 5X: lambda_X(n) = - sum_{e in {1,2,3,4}, e | n, n/e > X} mu(n/e).  Only X <= XMAX here (a falsification test, not a proof).
# uv run --no-project --with numpy python3 Riemann/Deep/R2AAudit/b3check.py
import numpy as np, time, resource
XMAX = 4*10**6; N = 5*XMAX
t0 = time.time()
mu = np.ones(N + 1, dtype=np.int8); isp = np.ones(N + 1, dtype=bool); isp[:2] = False
for i in range(2, int(N**0.5) + 1):
    if isp[i]: isp[i*i::i] = False
for p in np.nonzero(isp)[0]:
    p = int(p); mu[p::p] *= -1
    if p*p <= N: mu[p*p::p*p] = 0
del isp
mu[0] = 0
print('mu sieve to %d: %.1f s' % (N, time.time() - t0), flush=True)
for X in [10**4, 3*10**4, 10**5, 3*10**5, 10**6, 2*10**6, 4*10**6]:
    S = 0.0
    for a in range(X + 1, 5*X, 500000):                 # chunks of 5e5 (memory cap: every array <= 4 MB)
        n = np.arange(a, min(a + 500000, 5*X), dtype=np.int64)
        lam = -mu[n].astype(np.int64)
        for e in (2, 3, 4):
            m = (n % e == 0) & (n//e > X)
            lam[m] -= mu[n[m]//e]
        S += float(np.sum((lam*lam)/(n.astype(np.float64)**2)))
    print('X = %8d: X * sum_{X<n<5X} lambda_X(n)^2/n^2 = %.6f   (b3 = 0.605)' % (X, X*S), flush=True)
print('peak RSS %.0f MB, %.0f s' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20, time.time() - t0))
