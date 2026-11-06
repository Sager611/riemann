#!/usr/bin/env python3
"""d_chaingas.py -- explore15, transfer B (N7, N10, N12): hard-core gases at negative fugacity over the integers.

B1 (local lemma over primes).  Events A_p = {p | n} for n <= x: A_p, A_q exclusive when pq > x, so the
    dependency graph G_x on primes p <= x has an edge p ~ q iff pq > x (a threshold graph).  Shearer/Scott-Sokal:
    w = (p^-sigma) is in the Shearer region iff Z_{G_x}(-t w) > 0 for all t in [0, 1].  We compute the largest
    sigma where this fails, sigma_s(x), and the clique bound sigma_cl(x): sum_{sqrt x < q <= x} q^-sigma = 1.
B2 (chain gas, Knill/Bjorner).  H_n = incomparability graph of the squarefree integers in [2, n] under
    divisibility; independent sets = divisor chains, so Z_{H_n}(lam) = sum_{m <= n} mu^2(m) F_{omega(m)}(lam),
    F_j(lam) = sum_k k! S(j,k) lam^k (Fubini).  Philip Hall: Z_{H_n}(-1) = M(n).  We check the identity, locate
    the Shearer point lam_s(n) (first zero of Z(-lam)), and the zeros of Z_{H_n}.
Run: uv run --no-project --with numpy --with mpmath python3 -u d_chaingas.py   (about 1 min, < 300 MB)
"""
import math, time
import numpy as np
import mpmath as mp

def sieve(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]

def stirling2(J):
    S = [[0] * (J + 1) for _ in range(J + 1)]; S[0][0] = 1
    for j in range(1, J + 1):
        for k in range(1, j + 1): S[j][k] = k * S[j - 1][k] + S[j - 1][k - 1]
    return S

def main():
    t0 = time.time()
    N = 10 ** 7
    P = sieve(N)
    print(f'{len(P)} primes <= {N:.0e} ({time.time() - t0:.1f} s)')
    # ---- B1: Shearer threshold of the pair-knapsack dependency graph
    print('B1  x | sigma_cl(x) (clique of primes in (sqrt x, x]) | sigma_s(x) (Shearer, full graph G_x)')
    for e in range(2, 8):
        x = 10 ** e
        pr = P[P <= x].astype(float)
        big = pr[pr > math.sqrt(x)]
        f = lambda s: np.sum(big ** (-s)) - 1
        lo, hi = 0.2, 3.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if f(mid) > 0: lo = mid
            else: hi = mid
        scl = (lo + hi) / 2
        # threshold graph: independent sets = {} u {p_k} u subsets of primes <= min(p_k - 1, x/p_k)
        logp = np.log(pr)
        def Zt(s, t):
            w = t * pr ** (-s)
            c = np.cumsum(np.log1p(-w))                       # prod_{q <= p_i} (1 - t q^-s)
            lim = np.minimum(pr - 1, x / pr)
            idx = np.searchsorted(pr, lim, side='right') - 1  # last prime <= lim
            pre = np.where(idx >= 0, np.exp(c[np.maximum(idx, 0)]), 1.0)
            return 1 - np.sum(w * pre)
        def in_shearer(s):
            return all(Zt(s, t) > 0 for t in np.linspace(0.02, 1, 50))
        lo, hi = 0.2, 3.0
        for _ in range(40):
            mid = (lo + hi) / 2
            if in_shearer(mid): hi = mid
            else: lo = mid
        print(f'   1e{e} | {scl:.5f} | {hi:.5f}')
    print('    (both -> 1 as x -> oo: the local lemma over primes reaches only the abscissa of absolute convergence;'
          ' over F_q(t) the same argument gives |u| < 1/q, while RH lives on |u| = q^-1/2)')
    # ---- B2: chain gas
    t1 = time.time()
    omega = np.zeros(N + 1, dtype=np.int8)
    sqf = np.ones(N + 1, dtype=bool); sqf[0] = False
    for p in P:
        omega[p::p] += 1
        if p * p <= N: sqf[p * p::p * p] = False
    mu = np.where(sqf, np.where(omega % 2 == 0, 1, -1), 0).astype(np.int64)
    Mcum = np.cumsum(mu)
    print(f'B2  sieve to {N:.0e} ({time.time() - t1:.1f} s)')
    J = 9; S = stirling2(J)
    print('    n | M(n) | Z_{H_n}(-1) | #vertices V | lam_s (first zero of Z(-lam)) | lam_s*V | zeros of Z_{H_n}(lam)')
    for e in range(2, 8):
        n = 10 ** e
        cnt = np.bincount(omega[1:n + 1][sqf[1:n + 1]], minlength=J + 1)   # N_j, includes m = 1 (j = 0)
        c = [sum(int(cnt[j]) * math.factorial(k) * S[j][k] for j in range(J + 1)) for k in range(J + 1)]
        while c[-1] == 0: c.pop()
        Zm1 = sum(ck * (-1) ** k for k, ck in enumerate(c))
        mp.mp.dps = 50
        roots = mp.polyroots([mp.mpf(v) for v in reversed(c)], maxsteps=500, extraprec=300)
        real = [r for r in roots if abs(mp.im(r)) < 1e-30]
        cplx = [r for r in roots if abs(mp.im(r)) >= 1e-30]
        negreal = sorted([float(mp.re(r)) for r in real])
        lam_s = -max(negreal) if negreal else float('nan')
        V = int(cnt[1:].sum())
        zs = ', '.join(f'{v:.4f}' for v in negreal) + (f'; {len(cplx)} complex: ' + ', '.join(f'{complex(r):.3f}' for r in cplx) if cplx else '')
        print(f'   1e{e} | {int(Mcum[n])} | {Zm1} | {V} | {lam_s:.3e} | {lam_s * V:.3f} | {zs}')
    print(f'done ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
