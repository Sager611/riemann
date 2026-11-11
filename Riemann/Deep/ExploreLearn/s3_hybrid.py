# explore31 s3: where does the continuous relaxation gain?  Designs (indices k <-> dilation Q/k):
#   natural  {Q j : j <= N}                       (= Baez-Duarte 1/j)
#   head(K0) natural U {k : Q <= k <= Q K0}        (all theta = Q/k in [1/K0, 1] added)
#   full     {k : Q <= k <= Q N}                   (q-design, -> continuous relaxation as Q -> inf)
# Gram <rho_{Q/j}, rho_{Q/k}> = Q G_{jk}, G = Vasyunin (on-the-fly, cached V sums); target Q t_k, t_k = (log(k/Q)+1-gamma)/k.
# usage: uv run --no-project --with numpy --with scipy python3 s3_hybrid.py
import sys, math, time
import numpy as np
from scipy.linalg import cho_factor, cho_solve
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
LG = math.log(2 * math.pi) - EG
Vc = {}
def V(a, k):  # sum_{m=1}^{k-1} {a m/k} cot(pi m/k), a coprime to k
    if k == 1:
        return 0.0
    a %= k
    key = (a, k)
    if key in Vc:
        return Vc[key]
    h = (k - 1) // 2
    m = np.arange(1, h + 1)
    v = float(((2 * ((a * m) % k) - k) @ (1.0 / np.tan(np.pi * m / k))) / k) if h else 0.0
    if k <= 1500:
        Vc[key] = v
    return v
def vas(j, k):
    d = math.gcd(j, k); h, kk = j // d, k // d
    if h > kk:
        h, kk = kk, h
    return (LG * (h + kk) + (kk - h) * math.log(h / kk) - math.pi * (V(h, kk) + V(kk % h, h))) / (2.0 * h * kk) / d
def d2(idx, Q):
    idx = np.array(sorted(set(idx)))
    n = len(idx)
    Gs = np.empty((n, n))
    for a in range(n):
        for b in range(a, n):
            Gs[a, b] = Gs[b, a] = vas(int(idx[a]), int(idx[b]))
    t = (np.log(idx / Q) + 1 - EG) / idx
    cf = cho_factor(Gs, lower=True, check_finite=False)
    return 1 - Q * float(t @ cho_solve(cf, t, check_finite=False)), n
t0 = time.time()
RUNS = [(10, 64, [2, 3, 5, 10]), (40, 32, [2, 3, 5, 10, 40]), (160, 16, [2, 3, 5, 10, 20]), (640, 8, [2, 3, 5, 10]), (640, 32, [2, 3])]
if len(sys.argv) > 1 and sys.argv[1] == 'qext':  # head theta in [1/5,1], Q -> inf
    RUNS = [(N, Q, [5]) for N in (40, 160, 640) for Q in (8, 16, 32, 64)]
for N, Q, K0s in (RUNS if __name__ == '__main__' else []):
    L = math.log(N)
    nat = [Q * j for j in range(1, N + 1)]
    v0, n0 = d2(nat, Q)
    print(f'N={N} Q={Q}: natural d2L={v0*L:.6f} C2={L*(v0*L-C):+.5f} (n={n0})', flush=True)
    for K0 in K0s:
        v, n = d2(nat + list(range(Q, Q * K0 + 1)), Q)
        print(f'   head theta in [1/{K0},1]: d2L={v*L:.6f} C2={L*(v*L-C):+.5f} gain L^2(d2_nat-d2)={L*L*(v0-v):.5f} (n={n})  t={time.time()-t0:.0f}s', flush=True)
        check(250, f'N={N} K0={K0}')
if __name__ == '__main__':
    print('footprint/resident MB', mem(), 'cache', len(Vc))
