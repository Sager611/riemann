# T4 Pavlov transfer check (log|F_B| including sin(pi z); first run divided sin out and kept the lattice ripple): for a block B of unfolded zeta zeros (block-completed by Z outside), the generating function
# F_B(z) = sin(pi z) prod_{j in B} (z - lam_j)/(z - j) satisfies on Im z = y:
#   log|F_B(lam(t)+iy)/sin(pi(lam(t)+iy))| = log|zeta(1/2 + y/L(t) + it)| + smooth(t) + (outside-block correction).
import numpy as np, mpmath as mp
from common import zeta_unfolded, theta, dtheta, PI
mp.mp.dps = 15
g, lam = zeta_unfolded()
idx = np.arange(1, len(lam) + 1)
def lhs(t, y, lo, hi):
    z = theta(t)/PI + 1.5 + 1j*y
    B = slice(lo, hi)
    return np.log(np.abs(np.sin(PI*z))) + np.sum(np.log(np.abs(z - lam[B]))) - np.sum(np.log(np.abs(z - idx[B])))
for c, y in [(2000, 1.0), (8000, 1.0), (18000, 1.0), (8000, 0.5), (8000, 0.25), (18000, 0.25), (18000, 0.1)]:
    lo, hi = c - 300, c + 300
    ts = np.linspace(g[c - 20], g[c + 20], 81)
    A = np.array([lhs(t, y, lo, hi) for t in ts])
    Lt = dtheta(ts)/PI
    Bv = np.array([float(mp.log(abs(mp.zeta(mp.mpc(0.5 + y/l, t))))) for t, l in zip(ts, Lt)])
    X = np.vstack([Bv, np.ones_like(ts), ts - ts.mean()]).T
    coef, res, *_ = np.linalg.lstsq(X, A, rcond=None)
    r = A - X @ coef
    print(f'height {g[c]:8.1f} y={y:4.2f} sigma=1/2+{y/Lt.mean():.3f}: slope {coef[0]:.4f}  rms resid {r.std():.4f}  sd(log|zeta|) {Bv.std():.4f}  corr {np.corrcoef(A, Bv)[0,1]:.5f}', flush=True)
