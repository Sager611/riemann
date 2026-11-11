# nb2nd t5: (a) large-n shape of the N-free head optimizer beta*_K (K = 160, 320, 640); (b) Sigma_K for all K, fits of the convergence.
import math, numpy as np, mpmath as mp
O = '/Users/adriansager/playground/riemann/.scratch/nb2nd/'
Z = np.load(O + 'G640.npz'); G = Z['G']; b = Z['b']; mu = Z['mu'].astype(float); NMAX = len(b)
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
mp.mp.dps = 20
a1 = float(mp.zeta(0, derivative=2) / mp.zeta(0) - (mp.zeta(0, derivative=1) / mp.zeta(0)) ** 2)
z = [0.0] + [float(mp.zeta(2 * k + 1)) for k in range(1, 200)]
def ell(n):
    s0 = -(math.log(2 * math.pi) * (1 + math.log(n) - EG) + a1) / n
    if n == 1: tail = sum((z[k] - 1) / (2 * k * (2 * k + 1)) for k in range(1, 199)) + (1 - math.log(2))
    else:
        x = 1.0 / n; tail = sum(z[k] * x ** (2 * k + 1) / (2 * k * (2 * k + 1)) for k in range(1, 60))
    return s0 - tail
ks = np.arange(1, NMAX + 1); lth = np.array([ell(int(n)) for n in ks])
sf = mu != 0
for K in [160, 320, 640]:
    bs = np.linalg.solve(G[:K, :K], lth[:K]); n = ks[:K]; s_ = sf[:K]
    prof = np.where(s_, bs / np.where(s_, mu[:K], 1), np.nan)
    edges = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 641]
    out = []
    for i in range(len(edges) - 1):
        m = s_ & (n >= edges[i]) & (n < edges[i + 1]) & (n <= K)
        if m.any(): out.append((edges[i], round(float(np.nanmean(prof[m])), 3), round(float(np.nanstd(prof[m])), 3)))
    tailm = s_ & (n > K / 2)
    print(f'K={K}: beta*/mu binned by n (n_lo, mean, sd):', out)
    print(f'   last decile (n > 0.9K) mean beta/mu = {np.nanmean(prof[s_ & (n > 0.9*K)]):.3f}; non-sqfree rms (n>K/2) {np.sqrt(np.mean(bs[(~s_) & (n > K/2)]**2)):.3f}; sum beta/n = {bs @ (1/n):.4f}')
# Sigma_K for all K via nested Cholesky of G with rhs l
Lc = np.linalg.cholesky(G); yv = np.linalg.solve(Lc, lth); Sig = np.cumsum(yv * yv)
np.save(O + 't5_sigma.npy', Sig)
Ks = np.array([20, 40, 80, 160, 320, 640])
print('Sigma_K:', [(int(k), round(float(Sig[k - 1]), 6)) for k in Ks])
# fits on K in [40, 640]: (i) S_inf - A/log K, (ii) S_inf - A K^-th (th fitted), (iii) a + b log log K (divergent)
KK = np.arange(40, 641); S = Sig[KK - 1]
def lsq(X): coef, res, *_ = np.linalg.lstsq(X, S, rcond=None); return coef, float(np.sqrt(np.mean((X @ coef - S) ** 2)))
c1, e1 = lsq(np.vstack([np.ones_like(KK, dtype=float), 1 / np.log(KK)]).T)
c3, e3 = lsq(np.vstack([np.ones_like(KK, dtype=float), np.log(np.log(KK))]).T)
best = None
for th in np.linspace(0.1, 1.5, 57):
    c2, e2 = lsq(np.vstack([np.ones_like(KK, dtype=float), KK ** (-th)]).T)
    if best is None or e2 < best[2]: best = (th, c2, e2)
print(f'fit S_inf - A/log K: S_inf = {c1[0]:.5f}, A = {-c1[1]:.4f}, rms {e1:.2e}  -> J - S_inf = {J - c1[0]:.5f}')
print(f'fit S_inf - A K^-th: th = {best[0]:.3f}, S_inf = {best[1][0]:.5f}, rms {best[2]:.2e} -> J - S_inf = {J - best[1][0]:.5f}')
print(f'fit a + b loglog K (divergent): b = {c3[1]:.4f}, rms {e3:.2e}')
print(f'J - Sigma_640 = {J - Sig[-1]:.6f}')
