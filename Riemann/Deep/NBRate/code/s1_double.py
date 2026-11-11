# nbrate stage 1 (double, exploratory): Vasyunin Gram G (L^2(0,inf)), d_N^2 for all N (unconstrained V1 and
# constrained V2 = 2611.03924's d_n(1/2)^2), eigenvalues, BCF choice c_k = mu(k)(1 - log k/log N).
# usage: uv run --no-project --with numpy python3 s1_double.py NMAX
import sys, math, time
import numpy as np
from math import gcd
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1280
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
OUT = '/Users/adriansager/playground/riemann/.scratch/nbrate/'
# V(a,k) = sum_{m=1}^{k-1} {ma/k} cot(pi m/k) = (1/k) sum_{m<=(k-1)/2} (2(ma mod k) - k) cot(pi m/k), gcd(a,k)=1
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
def g(h, k):  # coprime
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
print(f'Gram built {time.time()-t0:.1f}s; G11={G[0,0]:.15f} G12={G[0,1]:.15f} G23={G[1,2]:.15f}', flush=True)
np.save(OUT + 'G_double.npy', G); np.save(OUT + 'b_double.npy', b)
# V1: unconstrained, all N via nested Cholesky
L = np.linalg.cholesky(G)
y = np.linalg.solve(L, b)  # general solve ok (triangular)
d1 = 1 - np.cumsum(y * y)  # d1[N-1] = d_N^2
# V2: constrained basis phi_k = rho_{1/k} - rho_1/k, k = 2..n ; A = T^T G T
A = G[1:, 1:] - np.outer(G[1:, 0], 1 / ks[1:]) - np.outer(1 / ks[1:], G[0, 1:]) + G[0, 0] * np.outer(1 / ks[1:], 1 / ks[1:])
b2 = np.log(ks[1:]) / ks[1:]
L2 = np.linalg.cholesky(A)
y2 = np.linalg.solve(L2, b2)
d2 = np.concatenate([[1.0], 1 - np.cumsum(y2 * y2)])  # d2[N-1], N>=1 (N=1: span empty -> 1)
# check 2611.03924's series for A_{jk} = (1/jk) sum_m (m mod j)(m mod k)/(m(m+1)) at small j,k (digamma over a period)
from math import lcm
def A_series(j, k, M=None):
    Lp = lcm(j, k)
    # sum over m>=1 of f(m)/(m(m+1)), f periodic with period Lp: sum_r f(r) sum_q 1/((qL+r)(qL+r+1))
    import mpmath as mp
    s = mp.mpf(0)
    for r in range(1, Lp + 1):
        f = (r % j) * (r % k)
        if f:
            s += f * (mp.digamma(mp.mpf(r + 1) / Lp) - mp.digamma(mp.mpf(r) / Lp)) / Lp
    return float(s / (j * k))
try:
    for (j, k) in [(2, 2), (2, 3), (3, 5), (4, 6), (5, 7), (6, 10)]:
        print(f'  A[{j},{k}] Vasyunin-based {A[j-2,k-2]:.15f}  series {A_series(j,k):.15f}', flush=True)
except Exception as e:
    print('series check skipped', e)
mu = np.ones(NMAX + 1, dtype=np.int64); mu[0] = 0
for p in range(2, NMAX + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu[p::p] *= -1
        mu[p * p::p * p] = 0
bcf = np.zeros(NMAX)
for N in range(2, NMAX + 1):
    c = mu[1:N + 1] * (1 - np.log(ks[:N]) / math.log(N))
    bcf[N - 1] = 1 - 2 * b[:N] @ c + c @ G[:N, :N] @ c
np.save(OUT + 'd1_double.npy', d1); np.save(OUT + 'd2_double.npy', d2); np.save(OUT + 'bcf_double.npy', bcf)
print(f'C = {C:.12f}')
print('N  d1=dN^2(unconstr)  d1*logN  d2(constr,2611)  d2*logN  BCF  BCF*logN')
for N in [10, 20, 40, 50, 80, 100, 160, 200, 320, 640, 1280]:
    if N <= NMAX:
        lN = math.log(N)
        print(f'{N:5d} {d1[N-1]:.10f} {d1[N-1]*lN:.6f}  {d2[N-1]:.10f} {d2[N-1]*lN:.6f}  {bcf[N-1]:.8f} {bcf[N-1]*lN:.6f}', flush=True)
lN = np.log(ks.astype(float)); lN[0] = np.nan
for name, d in (('V1', d1), ('V2', d2)):
    below = np.where(d[1:] * lN[1:] < C)[0] + 2
    print(name, 'first N with d^2 log N < C:', below[:12], ' count below:', len(below), 'of', NMAX - 1)
    # margins near crossing
print('eigen:')
for N in [40, 80, 160, 320, 640, 1280]:
    if N <= NMAX:
        ev = np.linalg.eigvalsh(G[:N, :N])
        ev2 = np.linalg.eigvalsh(A[:N - 1, :N - 1])
        print(f'  N={N}: G lmin={ev[0]:.4e} lmax={ev[-1]:.4e} cond={ev[-1]/ev[0]:.3e} | A lmin={ev2[0]:.4e} cond={ev2[-1]/ev2[0]:.3e}', flush=True)
print(f'done {time.time()-t0:.1f}s')
