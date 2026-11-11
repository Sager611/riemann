# explore31 s8 (control): is E_LC implemented right?  Gaussian sketches of the same pool (features psi_i = sum_j S_ij phi_j,
# S_ij iid N(0,1)) are the setting where the deterministic equivalent is accurate; compare with E_LC, random subsets, natural.
# theta = 1 pinned (Schur complement as in s2b).  usage: uv run --no-project --with numpy --with scipy python3 s8_sketch.py MMAX REPS
import sys, math, time
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
OUT = '/Users/adriansager/playground/riemann/.scratch/explore31/'
EG = 0.57721566490153286061
G = np.load(OUT + 'G2880.npy', mmap_mode='r')
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 50
rng = np.random.default_rng(8)
def dist2(Gs, b, y2):
    cf = cho_factor(Gs, lower=True, check_finite=False)
    return y2 - float(b @ cho_solve(cf, b, check_finite=False))
t0 = time.time()
print('N Q M | natural | E_LC | Gaussian sketch mean (sd) | random subset mean | sketch/E_LC')
for N in [3, 5, 10, 20, 40, 80, 160]:
    Q = max(1, (MMAX - 1) // (N - 1))
    ks = np.arange(Q, Q * N + 1)
    Gam = Q * np.array(G[Q - 1:Q * N, Q - 1:Q * N])
    beta = Q * (np.log(ks / Q) + 1 - EG) / ks
    g00 = Gam[0, 0]; g10 = Gam[1:, 0].copy(); b0 = beta[0]
    Gp = Gam[1:, 1:] - np.outer(g10, g10) / g00; bp = beta[1:] - g10 * b0 / g00; y2 = 1 - b0 * b0 / g00
    del Gam
    lam = eigh(Gp, eigvals_only=True, check_finite=False)
    lam2, U = eigh(Gp, check_finite=False); c = (U.T @ bp) ** 2; del U
    lam = np.maximum(lam2, 1e-300)
    n = N - 1
    lo, hi = -60.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if float(np.sum(lam / (lam + math.exp(mid)))) > n else (lo, mid)
    ELC = y2 - float(np.sum(c / (lam + math.exp(0.5 * (lo + hi)))))
    nat_idx = Q * np.arange(2, N + 1) - Q - 1
    nat = dist2(Gp[np.ix_(nat_idx, nat_idx)], bp[nat_idx], y2)
    sk, rs = [], []
    for r in range(REPS):
        S = rng.standard_normal((n, len(bp)))
        SG = S @ Gp
        sk.append(dist2(SG @ S.T, S @ bp, y2))
        idx = np.sort(rng.choice(len(bp), size=n, replace=False))
        rs.append(dist2(Gp[np.ix_(idx, idx)], bp[idx], y2))
    sk, rs = np.array(sk), np.array(rs)
    check(250, f'N={N}')
    print(f'{N:4d} {Q:4d} {len(ks):5d} | {nat:.6f} | {ELC:.6f} | {sk.mean():.6f} ({sk.std():.6f}) | {rs.mean():.6f} | {sk.mean()/ELC:.3f}  t={time.time()-t0:.0f}s', flush=True)
print('footprint/resident MB', mem())
