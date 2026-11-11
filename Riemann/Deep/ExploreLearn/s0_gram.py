# explore31 s0: double-precision Vasyunin Gram G_M = <rho(1/(hx)), rho(1/(kx))> in L^2(0,inf) (Baez-Duarte), h,k <= M.
# Copied from Riemann/Deep/NB2nd/t0_gram.py (itself from nbrate s1_double.py); G fill vectorised row-wise; output -> .scratch/explore31/.
# usage: uv run --no-project --with numpy python3 s0_gram.py M
import sys, math, time
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
t0 = time.time()
M = int(sys.argv[1]) if len(sys.argv) > 1 else 640
EG = 0.57721566490153286061
OUT = '/Users/adriansager/playground/riemann/.scratch/explore31/'
__import__('os').makedirs(OUT, exist_ok=True)
V = np.zeros((M + 1, M + 1))  # V[a,k] = sum_{m<k/2} (2{am/k}-1) cot(pi m/k) = sum_{m=1}^{k-1} {am/k} cot(pi m/k)
buf = np.empty(M * M // 2 + M, dtype=np.int64)  # reused: no allocator churn
fbuf = np.empty(M * M // 2 + M)
for k in range(2, M + 1):
    h = (k - 1) // 2
    if h == 0:
        continue
    m = np.arange(1, h + 1)
    a = np.array([x for x in range(1, k) if math.gcd(x, k) == 1])
    R = buf[:len(a) * h].reshape(len(a), h)
    np.multiply.outer(a, m, out=R); np.remainder(R, k, out=R)
    F = fbuf[:len(a) * h].reshape(len(a), h)
    np.multiply(R, 2.0, out=F); F -= k
    V[a, k] = (F @ (1.0 / np.tan(np.pi * m / k))) / k
    if k % 100 == 0:
        check(280, f'V k={k}')
del buf, fbuf, R, F
LG = math.log(2 * math.pi) - EG
G = np.zeros((M, M))
for hh in range(1, M + 1):
    kk = np.arange(hh, M + 1)
    d = np.gcd(hh, kk)
    h = hh // d; k = kk // d
    vhk = np.where(k > 1, V[h, k], 0.0)
    vkh = np.where(h > 1, V[k % np.maximum(h, 1), h], 0.0)
    row = (LG * (h + k) + (k - h) * np.log(h / k) - math.pi * (vhk + vkh)) / (2.0 * h * k) / d
    G[hh - 1, hh - 1:] = row
    G[hh - 1:, hh - 1] = row
del V
ks = np.arange(1, M + 1)
b = (np.log(ks) + 1 - EG) / ks
np.save(OUT + f'G{M}.npy', G)
f, r = mem()
print(f'G{M} built {time.time()-t0:.1f}s footprint {f:.0f}MB resident {r:.0f}MB; G11={G[0,0]:.15f} G12={G[0,1]:.15f} G23={G[1,2]:.15f}')
