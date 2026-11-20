# Checks for V14Twist: (1) k = thr < 1/2 along the chain (z < 7/6 - 3r/2); (2) finite zero-pad chain (step 1/4000, k = thr, wb = R) untwisted vs twisted;
# (3) break-even family saving gamma*(k) = (1-2k)/3 for a twisted-family estimate saving Z^(gamma tau) over the Cauchy count Z^tau.
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from v14_import import make
def thr(r, z): return (1.25 * r - 5 / 6) / (r - 0.5 - z)
G0, _ = make('free', 0.0)
import numpy as np
rs = np.linspace(0.7335, 7 / 9 - 1e-9, 20001)
zs = np.array([float(G0(r)) for r in rs]); ks = np.array([thr(r, z) for r, z in zip(rs, zs)])
print("(1) max over r in [0.7335, 7/9) of z(r) - (7/6 - 3r/2) = %.4e (negative <=> thr < 1/2); thr range [%.6f, %.6f]" % (max(zs - (7/6 - 1.5*rs)), ks.min(), ks.max()))
def chain(sig, step=1/4000, r_end=0.7335):
    R, z, out = 7 / 9, 0.0, {}
    while R > r_end:
        r = R - step; H = r / R; k = thr(R, z)
        e = (1 - H) / (6 + 2 * sig)            # 1 - H = 6e + 2 tau, tau = sig e
        z = H * z + e * (1 / k - 2); R = r
        for t in (0.7523, 0.7335):
            if t not in out and R <= t + 1e-12: out[t] = z
    return out
rho = 24 / 289; rp = rho / 3; sig = rp / (1 - rp)
for s in (0.0, sig, 3 * sig):
    o = chain(s); print("(2) finite chain step 1/4000 sig=%.6f: z(0.7523)=%.6e z(0.7335)=%.6e" % (s, o[0.7523], o[0.7335]))
g = (1 - 2 * ks) / 3
i = int(np.argmin(abs(rs - 0.752322))); print("(3) gamma* = (1-2k)/3: range [%.4f, %.4f]; at the contact crossing r=0.752322: k=%.6f gamma*=%.5f" % (g.min(), g.max(), ks[i], g[i]))
