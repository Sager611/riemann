# t5_posnb.py (explore22): positivity-constrained ('DRT-style', c_k >= 0) Baez-Duarte distance vs the unconstrained one.
# d_N^2 = min_c || chi_(0,1) - sum_{k<=N} c_k {1/(k t)} ||^2_{L^2(0,inf)} = 1 - 2 b.c + c.G c  (Vasyunin Gram, as explore13 nbgram.py).
# Elementary bound for c >= 0: on t > 1 the error is S/t, S = sum c_k/k; on (1/2,1) it is 1 + c_1 - S/t; hence
# d^2 >= min_{S, A>=1} [S^2 + A^2/2 - 2 ln2 A S + S^2] = (1 - ln^2 2)/2 = 0.259775 for every N.
import sys, math
import numpy as np
from math import gcd
EG = 0.5772156649015329
def V(h, k):
    if k == 1: return 0.0
    m = np.arange(1, k)
    return float(np.sum(((m * h) % k) / k / np.tan(np.pi * m / k)))
def gz_coprime(h, k):
    return ((math.log(2 * math.pi) - EG) * (h + k) / (2 * h * k) + (k - h) / (2 * h * k) * math.log(h / k)
            - math.pi / (2 * h * k) * (V(h, k) + V(k, h)))
def gram_zeta(N):
    G = np.zeros((N, N)); cache = {}
    for h in range(1, N + 1):
        for k in range(h, N + 1):
            g = gcd(h, k); a, b = h // g, k // g
            if (a, b) not in cache: cache[(a, b)] = gz_coprime(a, b)
            G[h - 1, k - 1] = G[k - 1, h - 1] = cache[(a, b)] / g
    b = np.array([(math.log(k) + 1 - EG) / k for k in range(1, N + 1)])
    return G, b
def nnls_cd(G, b, sweeps=20000, tol=1e-15, lo=None, hi=None):
    # coordinate descent for min c.Gc - 2b.c with box constraints lo <= c <= hi (default c >= 0)
    n = len(b); c = np.zeros(n)
    lo = np.zeros(n) if lo is None else lo; hi = np.full(n, np.inf) if hi is None else hi
    Gc = np.zeros(n)
    for it in range(sweeps):
        delta = 0.0
        for i in range(n):
            new = min(hi[i], max(lo[i], c[i] + (b[i] - Gc[i]) / G[i, i]))
            d = new - c[i]
            if d != 0.0:
                Gc += d * G[:, i]; c[i] = new; delta = max(delta, abs(d))
        if delta < tol: break
    return c, 1 - 2 * b @ c + c @ G @ c, it + 1
# DH dictionary (explore13): F(x) = A(x) = sum_{n<=x} a_n, A = 0 on [0,1): no 1/t tail (no pole), so the bound above does not apply.
kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
S5 = [0.0, 1.0, 1.0 + kap, 1.0, 0.0]
def gdh_coprime(h, k):
    from scipy.special import digamma
    P = 5 * h * k
    br = np.union1d(np.arange(0, P + 1, h), np.arange(0, P + 1, k)).astype(np.float64)
    a, bb = br[:-1], br[1:]
    c = np.array(S5)[(a // h).astype(np.int64) % 5] * np.array(S5)[(a // k).astype(np.int64) % 5]
    m = c != 0
    return float(np.sum(c[m] / P * (digamma(bb[m] / P) - digamma(a[m] / P))))
def gram_dh(N):
    from scipy.special import digamma
    G = np.zeros((N, N)); cache = {}
    for h in range(1, N + 1):
        for k in range(h, N + 1):
            g = gcd(h, k); a, b2 = h // g, k // g
            if (a, b2) not in cache: cache[(a, b2)] = gdh_coprime(a, b2)
            G[h - 1, k - 1] = G[k - 1, h - 1] = cache[(a, b2)] / g
    beta = sum(S5[r] * (digamma((r + 1) / 5) - digamma(r / 5)) / 5 for r in range(1, 5))
    return G, np.array([beta / k for k in range(1, N + 1)])
def d2_free(G, b, N):
    L = np.linalg.cholesky(G[:N, :N]); y = np.linalg.solve(L, b[:N]); return 1.0 - float(y @ y)
bound = (1 - math.log(2) ** 2) / 2
print('elementary lower bound (1 - ln^2 2)/2 =', round(bound, 6))
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 160
G, b = gram_zeta(Nmax)
for name, (G, b) in (('zeta', (G, b)), ('dh', gram_dh(min(Nmax, 80)))):
    for N in (1, 2, 5, 10, 20, 40, 80, 160):
        if N > len(b): break
        c, d2, its = nnls_cd(G[:N, :N], b[:N])
        supp = [k + 1 for k in range(N) if c[k] > 1e-12]
        print(name, 'N', N, 'd+^2 =', round(d2, 6), ' free d^2 =', round(d2_free(G, b, N), 6), ' support', supp[:14], 'sweeps', its, flush=True)

# Sign-pattern variants for zeta (which signs matter?):
#  (i) c_1 free, c_k >= 0 for k >= 2;  (ii) Mobius pattern: c_k * (-mu(k)) >= 0, c_k = 0 if mu(k) = 0;  (iii) c_1 <= 0, rest free.
def mobius(n):
    m, k, res = n, 2, 1
    while k * k <= m:
        if m % k == 0:
            m //= k
            if m % k == 0: return 0
            res = -res
        k += 1
    return -res if m > 1 else res
G, b = gram_zeta(Nmax)
for N in (10, 40, 160):
    if N > Nmax: break
    inf = np.inf
    lo1 = np.zeros(N); lo1[0] = -inf
    _, d_i, _ = nnls_cd(G[:N, :N], b[:N], lo=lo1, hi=np.full(N, inf))
    mu = np.array([mobius(k) for k in range(1, N + 1)], dtype=float)
    lo2 = np.where(mu < 0, 0.0, np.where(mu > 0, -inf, 0.0)); hi2 = np.where(mu < 0, inf, 0.0)
    _, d_ii, _ = nnls_cd(G[:N, :N], b[:N], lo=lo2, hi=hi2)
    lo3 = np.full(N, -inf); hi3 = np.full(N, inf); hi3[0] = 0.0
    _, d_iii, _ = nnls_cd(G[:N, :N], b[:N], lo=lo3, hi=hi3)
    print('zeta N', N, ' (i) c1 free, rest >= 0:', round(d_i, 6), ' (ii) signs of -mu(k):', round(d_ii, 6), ' (iii) c1 <= 0 only:', round(d_iii, 6), ' free:', round(d2_free(G, b, N), 6), flush=True)
