# nb2nd t0: double-precision Vasyunin Gram G (L^2(0,inf), V1) and b for N <= NMAX (copied from nbrate s1_double.py, trimmed).
# usage: uv run --no-project --with numpy python3 t0_gram.py NMAX   -> .scratch/nb2nd/G{NMAX}.npz
import sys, math, time
import numpy as np
from math import gcd
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 640
EG = 0.57721566490153286061
OUT = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
V = np.zeros((NMAX + 1, NMAX + 1))
for k in range(2, NMAX + 1):
    h = (k - 1) // 2
    if h == 0:
        continue
    m = np.arange(1, h + 1)
    a = np.array([x for x in range(1, k) if gcd(x, k) == 1])
    R = 2 * (np.outer(a, m) % k) - k
    V[a, k] = (R @ (1.0 / np.tan(np.pi * m / k))) / k
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
ks = np.arange(1, NMAX + 1)
b = (np.log(ks) + 1 - EG) / ks
mu = np.ones(NMAX + 1, dtype=np.int64); mu[0] = 0
for p in range(2, NMAX + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu[p::p] *= -1
        mu[p * p::p * p] = 0
np.savez(OUT + f'G{NMAX}.npz', G=G, b=b, mu=mu[1:])
print(f'G{NMAX} built {time.time()-t0:.1f}s; G11={G[0,0]:.15f} G12={G[0,1]:.15f} G23={G[1,2]:.15f}')
