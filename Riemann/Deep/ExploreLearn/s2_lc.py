# explore31 s2: the learning-curve ("omniscient risk", deterministic-equivalent) prediction vs actual designs.
# Pool = q-design {Q/k : Q <= k <= QN} (density of theta prop. 1/theta^2 on [1/N,1], the density of Baez-Duarte's 1/k).
# Pool Gram Gam = Q G[Q..QN], target beta = Q t, t_k = (log(k/Q)+1-gamma)/k.
# LC formula for n features sampled from the pool:  E_LC(n) = 1 - beta^T (Gam + eps)^{-1} beta,  tr(Gam (Gam+eps)^{-1}) = n.
# Compared with: natural d_N^2 (n = N, k = Q j), random n-subsets of the pool (Monte Carlo), full pool (~ continuous delta^2).
# usage: uv run --no-project --with numpy --with scipy python3 s2_lc.py MMAX REPS
import sys, math, time
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
OUT = '/Users/adriansager/playground/riemann/.scratch/explore31/'
__import__('os').makedirs(OUT, exist_ok=True)
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
G = np.load(OUT + 'G2880.npy', mmap_mode='r')
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 200
rng = np.random.default_rng(31)
def dist2(Gs, b):
    cf = cho_factor(Gs, lower=True, check_finite=False)
    return 1 - float(b @ cho_solve(cf, b, check_finite=False))
t0 = time.time()
print('N Q M | natural | E_LC(N) | random mean (sd, min) | pool(Q) | ratios natural/E_LC, random/E_LC')
for N in [2, 3, 5, 10, 20, 40, 80, 160, 320]:
    Q = max(1, (MMAX - 1) // (N - 1))
    while Q * N > G.shape[0]:
        Q -= 1
    ks = np.arange(Q, Q * N + 1)
    Gam = Q * np.array(G[Q - 1:Q * N, Q - 1:Q * N])
    beta = Q * (np.log(ks / Q) + 1 - EG) / ks
    lam, U = eigh(Gam, check_finite=False)
    c = (U.T @ beta) ** 2
    del U
    lam = np.maximum(lam, 1e-300)
    pool = 1 - float(np.sum(c / lam))
    lo, hi = -40.0, 10.0  # bisection on log eps for dof = N
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        dof = float(np.sum(lam / (lam + math.exp(mid))))
        lo, hi = (mid, hi) if dof > N else (lo, mid)
    eps = math.exp(0.5 * (lo + hi))
    ELC = 1 - float(np.sum(c / (lam + eps)))
    nat_idx = Q * np.arange(1, N + 1) - Q  # positions of k = Q j in ks
    nat = dist2(Gam[np.ix_(nat_idx, nat_idx)], beta[nat_idx])
    rs = []
    for r in range(REPS):
        idx = np.sort(rng.choice(len(ks), size=N, replace=False))
        rs.append(dist2(Gam[np.ix_(idx, idx)], beta[idx]))
    rs = np.array(rs)
    check(250, f'N={N}')
    L = math.log(N)
    print(f'{N:4d} {Q:4d} {len(ks):5d} | {nat:.6f} | {ELC:.6f} (eps={eps:.2e}) | {rs.mean():.6f} ({rs.std():.6f}, {rs.min():.6f}) | {pool:.6f} | {nat/ELC:.3f} {rs.mean()/ELC:.3f}  t={time.time()-t0:.0f}s', flush=True)
    print(f'      x logN: natural {nat*L:.5f}  E_LC {ELC*L:.5f}  random {rs.mean()*L:.5f}  pool {pool*L:.5f}   (C = {C:.5f})', flush=True)
print('footprint/resident MB', mem())
