# lp_ext.py: Krein/Bochner dual certificate test. Find an even extension E of W|(-L,L) on the gap L<=|x|<=P/2
# of the circle R/PZ such that every Fourier weight m(r_k)+E^(r_k) >= t (r_k = 2 pi k/P). t>0 => Q >= t ||f||^2.
# Usage: lp_ext.py c s basis(pt|hat) M B
import sys, math, numpy as np, time
from scipy.optimize import linprog
from weil_lat import *
c = float(sys.argv[1]); s = float(sys.argv[2]); basis = sys.argv[3]; M = int(sys.argv[4]); B = float(sys.argv[5])
L = math.log(c); at = atoms(c); P = 2 * L * (1 + s)
sig2 = 2 * sum(w for _, w in at)
# tail start: smallest r with hplus(r) - |T| - |pole| - sig2 - Ebound > 0 (monotone envelope)
def env(r, Eb):
    return hplus(r) - 2 * sum(math.exp(-(2*n+0.5)*L) * ((2*n+0.5) + r) / ((2*n+0.5)**2 + r*r) for n in range(40)) \
        - (2*math.exp(L/2)*2/ max(r - 0.5, 1e-9)) - sig2 - Eb
gap = P / 2 - L
if basis == 'pt':
    xs = np.linspace(L, P / 2, M)
    Eb = 2 * B
    def Ehat(r):
        return 2 * np.cos(np.outer(r, xs))
else:
    h = gap / (M + 1)
    xs = L + h * np.arange(1, M + 1)   # hats inside [L, P/2]; last one ends at P/2
    Eb = 2 * h * B
    def Ehat(r):
        r = np.asarray(r)
        sc = np.where(np.abs(r) < 1e-12, 1.0, np.sin(r * h / 2) / np.where(np.abs(r) < 1e-12, 1, r * h / 2)) ** 2
        return 2 * h * np.cos(np.outer(r, xs)) * sc[:, None]
rK = 10.0
while env(rK, Eb) <= 0.05:
    rK *= 1.05
K = int(math.ceil(rK * P / (2 * math.pi))) + 1
k = np.arange(0, K + 1); r = 2 * math.pi * k / P
m = m_sym(r, L, at)
A = Ehat(r)                       # (K+1) x M
nneg0 = int(np.sum(m < 0))
# variables: e+ (M), e- (M), t ; maximize t: minimize -t
nv = 2 * M + 1
cvec = np.zeros(nv); cvec[-1] = -1
# constraint: -(A e+ - A e-) + t <= m
Aub = np.hstack([-A, A, np.ones((K + 1, 1))]); bub = m.copy()
# sum(e+ + e-) <= B
Aub = np.vstack([Aub, np.concatenate([np.ones(2 * M), [0]])]); bub = np.concatenate([bub, [B]])
bounds = [(0, None)] * (2 * M) + [(None, 10)]
t0 = time.time()
res = linprog(cvec, A_ub=Aub, b_ub=bub, bounds=bounds, method='highs',
              options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
dt = time.time() - t0
if res.status != 0:
    print('LP status', res.status, res.message); sys.exit()
e = res.x[:M] - res.x[M:2 * M]; t = res.x[-1]
w = m + A @ e
act = int(np.sum(w - t < 1e-7 * max(1, abs(t))))
print(f'c={c:g} L={L:.4f} P={P:.4f} (s={s}) basis={basis} M={M} B={B}: K={K} nodes (r_K={rK:.1f}), '
      f'negative m(r_k) without E: {nneg0}; LP t* = {t:.4e}; min weight = {w.min():.4e}; '
      f'active nodes ~ {act}; sum|e| = {np.abs(e).sum():.3f}; LP {dt:.1f}s')
np.save(f'/Users/adriansager/playground/riemann/.scratch/explore8/ext_c{c:g}_s{s}_{basis}{M}.npy', np.concatenate([[c, s, P, t], xs, e]))
