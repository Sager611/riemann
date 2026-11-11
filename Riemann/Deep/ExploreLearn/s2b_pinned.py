# explore31 s2b: as s2, but theta = 1 (rho_1, the edge where chi jumps) is pinned in every design and the LC formula is
# applied to the remaining n = N-1 features after projecting out rho_1 (Schur complement): the fair version of the tool.
# usage: uv run --no-project --with numpy --with scipy python3 s2b_pinned.py MMAX REPS
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
def dist2(Gs, b, y2=1.0):
    cf = cho_factor(Gs, lower=True, check_finite=False)
    return y2 - float(b @ cho_solve(cf, b, check_finite=False))
t0 = time.time()
print('N Q M | natural | E_LC(N-1 | rho_1) | random+rho_1 mean (sd) | pool | natural/E_LC random/E_LC | dof-eps')
for N in [2, 3, 5, 10, 20, 40, 80, 160, 320]:
    Q = max(1, (MMAX - 1) // (N - 1))
    while Q * N > G.shape[0]:
        Q -= 1
    ks = np.arange(Q, Q * N + 1)
    Gam = Q * np.array(G[Q - 1:Q * N, Q - 1:Q * N])
    beta = Q * (np.log(ks / Q) + 1 - EG) / ks
    g00 = Gam[0, 0]; g10 = Gam[1:, 0].copy(); b0 = beta[0]
    Gp = Gam[1:, 1:] - np.outer(g10, g10) / g00
    bp = beta[1:] - g10 * b0 / g00
    y2 = 1 - b0 * b0 / g00
    lam, U = eigh(Gp, check_finite=False)
    c = (U.T @ bp) ** 2
    del U
    lam = np.maximum(lam, 1e-300)
    pool = y2 - float(np.sum(c / lam))
    n = N - 1
    lo, hi = -60.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        dof = float(np.sum(lam / (lam + math.exp(mid))))
        lo, hi = (mid, hi) if dof > n else (lo, mid)
    eps = math.exp(0.5 * (lo + hi))
    ELC = y2 - float(np.sum(c / (lam + eps)))
    nat_idx = Q * np.arange(2, N + 1) - Q - 1  # k = Q j, j >= 2, positions in the residual pool
    nat = dist2(Gp[np.ix_(nat_idx, nat_idx)], bp[nat_idx], y2) if n else y2
    rs = []
    for r in range(REPS):
        idx = np.sort(rng.choice(len(ks) - 1, size=n, replace=False))
        rs.append(dist2(Gp[np.ix_(idx, idx)], bp[idx], y2))
    rs = np.array(rs)
    check(250, f'N={N}')
    L = math.log(N)
    print(f'{N:4d} {Q:4d} {len(ks):5d} | {nat:.6f} | {ELC:.6f} | {rs.mean():.6f} ({rs.std():.6f}) | {pool:.6f} | {nat/ELC:.3f} {rs.mean()/ELC:.3f} | eps={eps:.2e}  t={time.time()-t0:.0f}s', flush=True)
    print(f'      x logN: natural {nat*L:.5f}  E_LC {ELC*L:.5f}  random {rs.mean()*L:.5f}  pool {pool*L:.5f}   (C = {C:.5f}); rho_1 alone {y2:.5f}', flush=True)
print('footprint/resident MB', mem())
