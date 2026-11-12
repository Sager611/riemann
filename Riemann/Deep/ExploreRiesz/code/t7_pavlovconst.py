# T7: long-window Pavlov constant of zeta zeros: [w]_{A2} for w = |zeta(sigma + it)|^2, sigma = 1/2 + y/L(T), vs zeta(2s)^2/zeta(4s)
# (Carlson: mean |zeta|^2 = zeta(2s), mean |zeta|^-2 = zeta(2s)/zeta(4s); unconditional for s > B = 0.87338, RH for s > 1/2).
import numpy as np, mpmath as mp
from common import dtheta, PI
mp.mp.dps = 15
def a2(ts, w, Kt):
    dt = ts[1] - ts[0]; m = max(1, int(round(Kt/dt)))
    c1 = np.concatenate([[0], np.cumsum(w)]); c2 = np.concatenate([[0], np.cumsum(1/w)])
    a = (c1[m:] - c1[:-m])/m; b = (c2[m:] - c2[:-m])/m; p = a*b
    return np.median(p), p.max()
for Tc, y in [(11000.0, 1.0), (11000.0, 0.5), (11000.0, 0.25)]:
    L = dtheta(Tc)/PI; s = 0.5 + y/L
    ts = np.arange(Tc - 220, Tc + 220, 0.05/L)
    w = np.array([float(abs(mp.zeta(mp.mpc(s, t))))**2 for t in ts])
    pred = float(mp.zeta(2*s)**2/mp.zeta(4*s))
    print(f'T~{Tc:.0f} y={y} sigma={s:.4f}: mean w = {w.mean():.4f} (zeta(2s) = {float(mp.zeta(2*s)):.4f}), mean 1/w = {(1/w).mean():.4f} (zeta(2s)/zeta(4s) = {float(mp.zeta(2*s)/mp.zeta(4*s)):.4f}), predicted A2 limit {pred:.4f}', flush=True)
    print('   [w]_A2 median/max over windows of K zeros: ' + '  '.join(f'K={K}: {a2(ts, w, K/L)[0]:.3f}/{a2(ts, w, K/L)[1]:.2f}' for K in [1, 4, 16, 64, 256]), flush=True)
