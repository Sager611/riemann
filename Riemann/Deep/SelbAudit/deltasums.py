#!/usr/bin/env python3
"""selbaudit check 1 (DERIVATION §10): signed local delta-sums, exact A_d ratios, exact Sigma_1 at finite xi.
uv run --no-project --with numpy python3 deltasums.py   (single thread, < 200 MB)"""
import math, sys, time
import numpy as np
T0 = time.time()
def tau_pk(a):  # tau_{-1/2}(p^a) = binom(a-3/2, a)
    r = 1.0
    for i in range(1, a + 1): r *= (i - 1.5)/i
    return r
def tau_half_pk(a):
    r = 1.0
    for i in range(1, a + 1): r *= (i - 0.5)/i
    return r
print("== (a) sigma_p(k,v) = sum_{a+b>=k} tau(p^a)tau(p^b) p^{-a(1+iv)-b}  (claim: -1/p (k=1), 0 (k>=2) at v=0) ==")
for p in [2, 3, 5, 7, 11, 13]:
    row = []
    for k in [1, 2, 3, 4]:
        for v in [0.0, 0.01, 0.1]:
            s = 0j
            for a in range(0, 80):
                for b in range(0, 80):
                    if a + b >= k: s += tau_pk(a)*tau_pk(b)*p**(-a - b)*complex(math.cos(a*v*math.log(p)), -math.sin(a*v*math.log(p)))
            row.append((k, v, s))
    absk = [sum(abs(tau_pk(a))*abs(tau_pk(b))*p**(-a - b) for a in range(80) for b in range(80) if a + b >= k) for k in [1, 2, 3, 4]]
    hk = [sum(tau_half_pk(a)*tau_half_pk(b)*p**(-a - b) for a in range(80) for b in range(80) if a + b >= k) for k in [1, 2, 3, 4]]
    print("p=%-3d " % p + " ".join("k%d,v%.2g:%+.6f%+.6fi" % (k, v, s.real, s.imag) for (k, v, s) in row if v == 0.0),
          "| -1/p=%.6f | |signed-sum|_abs(tau_-1/2) k=1..4: %s | tau_1/2-sum k=1..4: %s vs 1/phi(p^k): %s" % (-1/p,
          " ".join("%.6f" % x for x in absk), " ".join("%.6f" % x for x in hk), " ".join("%.6f" % (1/(p**(k - 1)*(p - 1))) for k in [1, 2, 3, 4])))
    print("      v-drift |sigma(k,v)-sigma(k,0)| at v=0.01,0.1: " + " ".join("k%d:%.2e,%.2e" % (k, abs(row[(k-1)*3+1][2]-row[(k-1)*3][2]), abs(row[(k-1)*3+2][2]-row[(k-1)*3][2])) for k in [1,2,3,4]))

def tau_sieve(n):
    t = np.ones(n + 1); spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, n + 1):
        if spf[i] == 0:
            spf[i::i][spf[i::i] == 0] = i
    t[0] = 0
    for m in range(2, n + 1):
        p = spf[m]; r = m; a = 0
        while r % p == 0: r //= p; a += 1
        t[m] = t[r]*tau_pk(a)
    return t, spf
def coprime_mask(n, primes):
    m = np.ones(n + 1, dtype=bool); m[0] = False
    for p in primes: m[::p] = False
    return m
print("== (b) exact A_d/A_1 at finite xi (q = 1 and q = 13), theta = 0.15; claim (signed, xi -> oo): (-1)^w(d) mu^2(d)/phi(d); absolute-sum bound: (d/phi(d))/phi(d)/Bhat(0)^2 ==")
TH = 0.15
def run_b(xi, q):
    n = xi; t, spf = tau_sieve(n)
    L = math.log(xi); Lt = TH*L
    idx = np.arange(n + 1, dtype=float); idx[0] = 1.0
    lg = np.log(idx)
    out = {}
    cache = {}
    def S(mask_primes, z):   # returns (S0(z), S1(z)) = sum_{n<z, cop} tau(n)/n, tau(n) log n / n
        key = tuple(sorted(mask_primes))
        if key not in cache:
            m = coprime_mask(n, key); w = np.where(m, t/idx, 0.0); w[0] = 0
            cache[key] = (np.cumsum(w), np.cumsum(w*lg))
        c0, c1 = cache[key]; k = int(math.ceil(z)) - 1   # n < z  <=> n <= ceil(z)-1
        k = min(max(k, 0), n)
        return c0[k], c1[k]
    def B(d_primes, qpr, delta, tdelta):
        y = xi/delta; yt = xi**TH/delta
        s0, s1 = S(d_primes + qpr, y); v = math.log(y)*s0 - s1
        if yt > 1:
            s0, s1 = S(d_primes + qpr, yt); v -= math.log(yt)*s0 - s1
        return tdelta*v/((1 - TH)*L)
    qpr = [] if q == 1 else [13]
    def A(d):
        f = {}; r = d; p = 2
        while r > 1:
            while r % p == 0: f[p] = f.get(p, 0) + 1; r //= p
            p += 1
        ps = sorted(f)
        # enumerate delta | d^infty, delta < xi
        dels = [(1, 1.0, {})]
        for p in ps:
            new = []
            for (dl, tl, ex) in dels:
                a = 0; v = dl
                while v < xi:
                    e2 = dict(ex); e2[p] = a; new.append((v, tl*tau_pk(a), e2)); v *= p; a += 1
            dels = new
        Bv = [(dl, B(ps, qpr, dl, tl), ex) for (dl, tl, ex) in dels]
        s = 0.0
        for (d1, b1, e1) in Bv:
            for (d2, b2, e2) in Bv:
                if all(e1[p] + e2[p] >= f[p] for p in ps): s += b1*b2/(d1*d2)
        return s
    A1 = A(1); phi = lambda d: sum(1 for k in range(1, d + 1) if math.gcd(k, d) == 1)
    Bh0 = (1 - math.sqrt(TH))/math.sqrt(1 - TH)
    res = []
    for d in [2, 3, 4, 5, 6, 8, 9, 10, 12, 30]:
        w = len([p for p in range(2, d + 1) if d % p == 0 and all(p % r for r in range(2, p))])
        sqf = all(d % (p*p) for p in range(2, d + 1))
        pred = ((-1)**w/phi(d)) if sqf else 0.0
        res.append((d, A(d)/A1, pred, (d/phi(d))/phi(d)/Bh0**2))
    return A1, res
for q in [1, 13]:
    for xi in [10**4, 10**5, 10**6]:
        A1, res = run_b(xi, q)
        print("q=%-2d xi=1e%d  A_1*log(xi)=%.5f | " % (q, round(math.log10(xi)), A1*math.log(xi)) + " ".join("d%d:%+.4f(pred %+.3f,abs %.2f)" % r for r in res) + "  (%.0fs)" % (time.time() - T0), flush=True)
