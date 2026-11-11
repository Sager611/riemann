# explore31 s7: does nb2nd's head functional explain the continuum head gain?
# l(x) = <zeta'/zeta, zeta x^{-s}>_mu for real x >= 1 (nb2nd eq. (2), residues at s = 0, -2k; x^{s-1} in place of n^{s-1}):
#   l(x) = -[log(2pi)(1 + log x - gamma) + a1]/x - sum_k (zeta(2k+1)-1) x^{-2k-1}/(2k(2k+1)) - [-(y/2)log(1-y^2) - atanh(y) + y], y = 1/x,
#   a1 = pi^2/12 - gamma^2 - 2 gamma_1.  Head design H (dilations 1/x), Sigma_H = l_H^T Gram_H^{-1} l_H, Gram <zeta x1^{-s}, zeta x2^{-s}> = Q G_{m1 m2}, x = m/Q.
# Prediction (nb2nd step 3 for a fixed head): L^2 (d_nat^2 - d_{nat U H}^2) -> Sigma_{nat U H} - Sigma_nat.
# usage: uv run --no-project --with mpmath --with numpy --with scipy python3 s7_headfun.py
import sys, math, time
import numpy as np
import mpmath as mp
from scipy.linalg import cho_factor, cho_solve
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from memcap import check, mem
from s3_hybrid import vas
EG = 0.57721566490153286061
a1 = float(mp.pi ** 2 / 12 - mp.euler ** 2 - 2 * mp.stieltjes(1))
L2P = math.log(2 * math.pi)
zk = [float(mp.zeta(2 * k + 1) - 1) for k in range(1, 60)]
def ell(x):
    y = 1.0 / x
    s = sum(z * y ** (2 * k + 1) / (2 * k * (2 * k + 1)) for k, z in enumerate(zk, start=1))
    br = (1 - math.log(2)) if x == 1 else (-(y / 2) * math.log(1 - y * y) - math.atanh(y) + y)
    return -(L2P * (1 + math.log(x) - EG) + a1) / x - s - br
print(f'l(1) = {ell(1.0):.10f} (nb2nd: Sigma_1 = l1^2/G11 = {ell(1.0)**2/1.260661401507813:.6f}, quoted 2.4419)')
def sigma(idx, Q):
    idx = sorted(set(idx)); n = len(idx)
    Gs = np.empty((n, n))
    for a in range(n):
        for b in range(a, n):
            Gs[a, b] = Gs[b, a] = Q * vas(idx[a], idx[b])
    l = np.array([ell(m / Q) for m in idx])
    cf = cho_factor(Gs, lower=True, check_finite=False)
    return float(l @ cho_solve(cf, l, check_finite=False)), n
t0 = time.time()
for K in (40, 80, 160):
    for Q in (8, 16, 32, 64):
        nat = [Q * j for j in range(1, K + 1)]
        s0, n0 = sigma(nat, Q)
        s5, n5 = sigma(nat + list(range(Q, 5 * Q + 1)), Q)
        s2, n2 = sigma(nat + list(range(Q, 2 * Q + 1)), Q)
        check(250, f'K={K} Q={Q}')
        print(f'K={K:4d} Q={Q:3d}: Sigma_nat={s0:.6f}  +[1/2,1]: {s2:.6f} (gain {s2-s0:.5f})  +[1/5,1]: {s5:.6f} (gain {s5-s0:.5f})  t={time.time()-t0:.0f}s', flush=True)
print('J = 3.5626634565; footprint/resident MB', mem())
