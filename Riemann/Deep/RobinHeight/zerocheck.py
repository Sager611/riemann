#!/usr/bin/env python3
# zerocheck.py (about 2-3 minutes; uv run --no-project --with numpy --with mpmath python3 zerocheck.py) -- numerical validation (not part of the proof) of Lemma 5's explicit formula, signs and constants:
#   -I(x) = M(x) + Z(x) + log(2 pi)/(x log x) + V(x),  Z(x) = sum_rho (1/rho) int_x^inf t^(rho-2) h(t) dt,
# with -I(x) computed from primes via Lemma 1/2 and M(x) from prime powers, Z from the first K zeros.
import numpy as np, mpmath as mp, math
mp.mp.dps = 30
N = 10**7
sv = np.ones(N + 1, dtype=bool); sv[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.float64)
logp = np.log(P)
cum_theta = np.cumsum(logp)
cum_logP = np.cumsum(-np.log1p(-1.0/P))
EG = float(mp.euler)
tailE1 = float(mp.e1(math.log(N)))
K = 300
zeros = [mp.zetazero(k) for k in range(1, K + 1)]
def Fz(rho, L):
    return rho*mp.e1((1 - rho)*L) + mp.exp(-(1 - rho)*L)/L
for x in [10**5, 10**6, 3*10**6]:
    j = np.searchsorted(P, x, side='right')          # primes <= x : P[:j]
    th = cum_theta[j-1]; lP = cum_logP[j-1]
    L = math.log(x)
    D = math.log(L) - math.log(math.log(th)) + (th - x)/(x*L)
    pg = P[j:]
    S2 = float(np.sum(-np.log1p(-1.0/pg) - 1.0/pg)) + tailE1/2/N*0   # tail beyond 1e7 ~ 1/(2*1e7*16): add below
    S2 += 1.0/(2*N*math.log(N))
    negI = (lP - EG - math.log(math.log(th))) - D + S2
    # M(x) = (psi-theta)(x)/(x log x) + sum_{k>=2} sum_{p^k > x} 1/(k p^k)
    psith = 0.0; Msum = 0.0
    for k in range(2, 60):
        r = x**(1.0/k)
        if r < 2: break
        jj = np.searchsorted(P, r, side='right')
        psith += cum_theta[jj-1] if jj > 0 else 0.0
        pk = P[jj:]
        if k == 2:
            Msum += float(np.sum(1.0/(2*pk*pk))) + 1.0/(2*N*math.log(N))   # tail p > 1e7
        else:
            Msum += float(np.sum(pk**(-float(k))))/k
    M = psith/(x*L) + Msum
    Ztrue = negI - M - math.log(2*math.pi)/(x*L)
    ZK = sum(2*mp.re(Fz(r, L)/r) for r in zeros)
    n = math.sqrt(x)*L
    print(f'x={x:.0e}: -I*n={negI*n:.6f}  M*n={M*n:.6f}  Ztrue*n={Ztrue*n:+.6f}  Z_K*n={float(ZK)*n:+.6f}  (K={K}, |diff|={abs(Ztrue-float(ZK))*n:.2e})')
