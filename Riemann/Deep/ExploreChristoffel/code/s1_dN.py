# explore38 s1: double-precision Vasyunin Gram (V1, L^2(0,inf)) and d_N^2 for every N <= NMAX via one Cholesky.
# Gram code copied from Riemann/Deep/NB2nd/t0_gram.py (itself from nbrate s1_double.py).
# usage: uv run --no-project --with numpy --with scipy python3 s1_dN.py NMAX
import sys, math, time
import numpy as np
import resource
def mem(tag): print(tag, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6, 'MB', flush=True)
from math import gcd
from scipy.linalg import cholesky, solve_triangular
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1280
EG = 0.57721566490153286061
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
V = np.zeros((NMAX + 1, NMAX + 1))
for k in range(2, NMAX + 1):
    h = (k - 1) // 2
    if h == 0:
        continue
    m = np.arange(1, h + 1)
    a = np.array([x for x in range(1, k) if gcd(x, k) == 1], dtype=np.int64)
    ct = 1.0 / np.tan(np.pi * m / k)
    for i0 in range(0, len(a), 32):  # chunked: the full outer product leaked ~390 MB at NMAX = 1280
        aa = a[i0:i0 + 32]
        R = (2 * ((aa[:, None] * m[None, :]) % k) - k).astype(np.float64)
        V[aa, k] = (R @ ct) / k
mem('after V')
LG = math.log(2 * math.pi) - EG
def g(h, k):  # coprime h, k
    vhk = V[h, k] if k > 1 else 0.0
    vkh = V[k % h, h] if h > 1 else 0.0
    return (LG * (h + k) + (k - h) * math.log(h / k) - math.pi * (vhk + vkh)) / (2 * h * k)
G = np.zeros((NMAX, NMAX))
for hh in range(1, NMAX + 1):
    for kk in range(hh, NMAX + 1):
        d = gcd(hh, kk)
        G[hh - 1, kk - 1] = G[kk - 1, hh - 1] = g(hh // d, kk // d) / d
del V
mem('after G')
ks = np.arange(1, NMAX + 1)
b = (np.log(ks) + 1 - EG) / ks
Lc = cholesky(G, lower=True, overwrite_a=True, check_finite=False)
y = solve_triangular(Lc, b, lower=True, check_finite=False)
d2 = 1.0 - np.cumsum(y * y)
mem('after chol')
np.savez(OUT + f'dN{NMAX}.npz', d2=d2)
C = 2 + EG - math.log(4 * math.pi)
for N in [2, 10, 40, 41, 50, 80, 160, 200, 320, 640, 1280]:
    if N <= NMAX:
        print(f'N={N:5d} d2={d2[N-1]:.15f} d2logN={d2[N-1]*math.log(N):.10f} minusC={d2[N-1]*math.log(N)-C:+.4e}')
print(f'done {time.time()-t0:.1f}s')
