# T6a: band-limited Wirtinger gain vs number m of consecutive maximal gaps (density 1 restored by equal small gaps).
# T6b: Pavlov A2 characteristic [w]_{A2} of w = |F_B(x+iy)|^2 over unfolded windows of length K: zeta, CUE, DH, Poisson.
import numpy as np, sys
from scipy.linalg import null_space
from common import zeta_unfolded, dh_unfolded, cue, PI
from t5_bandwirt import R_config
kap = 3.324; I = (-1.5*PI, PI); n = 40
print('T6a kappa=3.324, band [-1.5pi, pi], period n=40: m consecutive gaps kappa; R/(pi/kappa)^2')
for m in [1, 2, 3, 5, 8, 11]:
    rest = n - m*kap; k = n - m        # remaining n - m gaps share length n - m*kappa
    x = np.concatenate([kap*np.arange(m + 1), m*kap + rest*np.arange(1, k)/k])
    r = R_config(x, n, I, nq=121)
    print(f'  m={m:2d} small gap {rest/k:.3f}: ratio {r/(PI/kap)**2:.4f}', flush=True)
def a2(lam_blk, idx_blk, xs, y):
    z = xs + 1j*y
    lw = 2*np.log(np.abs(np.sin(PI*z)))
    for l, j in zip(lam_blk, idx_blk):
        lw += 2*(np.log(np.abs(z - l)) - np.log(np.abs(z - j)))
    return lw
def a2char(lw, dx, K):
    w = np.exp(lw - lw.mean()); c1 = np.concatenate([[0], np.cumsum(w)]); c2 = np.concatenate([[0], np.cumsum(1/w)])
    m = int(round(K/dx)); a = (c1[m:] - c1[:-m])/m; b = (c2[m:] - c2[:-m])/m
    p = a*b; return np.median(p), p.max()
rng = np.random.default_rng(7)
g, lam = zeta_unfolded()
sets = {}
c = 12000; sets['zeta(t~1e4)'] = (lam[c-800:c+800], np.arange(c-800, c+800) + 1.0, c - 256.0, c + 256.0)
cu = cue(1024, rng); sets['CUE'] = (cu + 100, np.arange(1024) + 100.5, 356.0, 868.0)
po = np.cumsum(rng.exponential(1, 1024)); po = po*1024/po[-1]; sets['Poisson'] = (po + 100, np.arange(1024) + 100.5, 356.0, 868.0)
pts, _ = dh_unfolded(True); sets['DH(T<420)'] = (pts, np.arange(1, len(pts) + 1) + 0.0, 60.0, 260.0)
print('\nT6b [w]_A2 median / max over windows of K unfolded units (line Im z = y)')
for y in [1.0, 0.5]:
    for name, (L_, J_, a, b) in sets.items():
        xs = np.arange(a, b, 0.05)
        lw = a2(L_, J_, xs, y)
        out = []
        for K in [1, 4, 16, 64, 128]:
            if K > (b - a)/2: continue
            md, mx = a2char(lw, 0.05, K); out.append(f'K={K}: {md:.3f}/{mx:.2f}')
        print(f'  y={y} {name:12s} ' + '  '.join(out), flush=True)
