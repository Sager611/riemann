# explore13 T2 test: Baez-Duarte (natural numbers, product order on Q+*, where Nehari/AAK fail) versus the
# half-space relaxation (all rational dilations theta = a/b in [1/N, 1], b <= B: the linearly ordered side where
# Helson-Lowdenslager/AAK hold).  d_BD(N) >= D(log N) (continuous Beurling problem with memory log N) <= d_Farey(N,B).
# Gram: <rho_alpha, rho_beta> = a c Vas(b c, a d) for alpha = a/b, beta = c/d (homogeneity), target theta(log(1/theta)+1-gamma).
# Usage: python3 farey.py N B dps
import sys
from math import gcd
import mpmath as mp
N = int(sys.argv[1]); B = int(sys.argv[2]); mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 40
cache = {}
def Vc(h, k):  # sum_{m=1}^{k-1} {m h/k} cot(pi m/k)
    return mp.fsum(mp.mpf((m * h) % k) / k * mp.cot(mp.pi * m / k) for m in range(1, k))
def vas(h, k):
    g = gcd(h, k); h //= g; k //= g
    if h > k: h, k = k, h
    if (h, k) not in cache:
        cache[(h, k)] = ((mp.log(2 * mp.pi) - mp.euler) * (h + k) / (2 * h * k) + mp.mpf(k - h) / (2 * h * k) * mp.log(mp.mpf(h) / k)
                         - mp.pi / (2 * h * k) * (Vc(h, k) + Vc(k, h)))
    return cache[(h, k)] / g
def dist2(dic):
    n = len(dic)
    G = mp.matrix(n, n); bb = mp.matrix(n, 1)
    for i, (a, b) in enumerate(dic):
        th = mp.mpf(a) / b
        bb[i] = th * (mp.log(1 / th) + 1 - mp.euler)
        for j in range(i, n):
            c, d = dic[j]
            G[i, j] = G[j, i] = a * c * vas(b * c, a * d)
    L = mp.cholesky(G)
    y = mp.lower_triangular_solve(L, bb) if hasattr(mp, 'lower_triangular_solve') else mp.lu_solve(L, bb)
    return 1 - mp.fsum(y[i] ** 2 for i in range(n)), n
bd = [(1, k) for k in range(1, N + 1)]
far = sorted({(a // gcd(a, b), b // gcd(a, b)) for b in range(1, B + 1) for a in range(1, b + 1) if a * N >= b}, key=lambda t: t[0] / t[1])
d_bd, n1 = dist2(bd)
d_far, n2 = dist2(far)
L = mp.log(N)
print(f"N={N} B={B} dps={mp.mp.dps}: BD size {n1} d^2={mp.nstr(d_bd, 10)} d^2 logN={mp.nstr(d_bd * L, 6)} | Farey size {n2} d^2={mp.nstr(d_far, 10)} d^2 logN={mp.nstr(d_far * L, 6)} | ratio {mp.nstr(d_far / d_bd, 5)}", flush=True)
