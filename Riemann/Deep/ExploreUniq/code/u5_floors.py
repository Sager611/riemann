# explore37 u5: floors eta_Mc(S)^2 at larger Mc, memory-lean (Gram accumulated over Gauss cylinders in chunks;
# no full V matrix). Same time-domain construction as u1_exact.py. usage: python3 u5_floors.py M
# (first version materialised V and fancy-index copies: peak 439 MB > 300 MB cap; replaced by this one.)
import sys, math, time, resource
import numpy as np, mpmath as mp
from scipy.linalg.blas import dsyrk
from scipy.linalg import cho_factor, cho_solve
t0 = time.time(); M = int(sys.argv[1]); mp.mp.dps = 40
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
R00 = np.zeros(M + 1); R01 = np.zeros(M + 1); R11 = np.zeros(M + 1)
for j in range(1, M + 1):
    a_ = mp.mpf(1) / (j + 1); b_ = mp.mpf(1) / j; Lj = mp.log(j)
    F1 = lambda x: x * (mp.log(x) + Lj) - x
    F2 = lambda x: (lambda lx: x * lx**2 - 2 * x * lx + 2 * x)(mp.log(x) + Lj)
    m0 = b_ - a_; m1 = F1(b_) - F1(a_); m2 = F2(b_) - F2(a_)
    r00 = mp.sqrt(m0); r01 = m1 / r00; r11 = mp.sqrt(m2 - r01**2)
    R00[j], R01[j], R11[j] = float(r00), float(r01), float(r11)
lgam = np.array([math.lgamma(n + 1) for n in range(M + 2)])
mu = np.ones(M + 1, dtype=np.int64); mu[0] = 0
for p in range(2, M + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu[p::p] *= -1; mu[p * p::p * p] = 0
def rows_block(j0, j1, ks):  # local orthonormal coords on cylinders j0..j1-1 for columns ks (array) and for e
    jj = np.arange(j0, j1); lj = np.log(jj.astype(float))
    n = jj[:, None] // ks[None, :]  # floor(j/k)
    alpha = n * (1 + np.log(ks.astype(float))[None, :] - lj[:, None]) + lgam[n]; beta = n.astype(float)
    B = np.empty((2 * len(jj), len(ks)))
    B[0::2] = R00[jj][:, None] * alpha + R01[jj][:, None] * beta; B[1::2] = R11[jj][:, None] * beta
    y = np.empty(2 * len(jj)); y[0::2] = R00[jj] * (1 - lj) + R01[jj]; y[1::2] = R11[jj]
    return B, y
def floors(S, Ms, chunk=50):
    ks = np.array(sorted(k for k in S if k <= max(Ms))); G = np.zeros((len(ks), len(ks)), order='F'); r = np.zeros(len(ks)); yy = 0.0
    out = []; j = 1
    for Mc in sorted(Ms):
        while j <= Mc:
            j1 = min(j + chunk, Mc + 1); B, y = rows_block(j, j1, ks)
            G = dsyrk(1.0, B, beta=1.0, c=G, trans=1, lower=0, overwrite_c=1); r += B.T @ y; yy += y @ y; j = j1
        m = int(np.searchsorted(ks, Mc, side='right'))
        cf = cho_factor(G[:m, :m], lower=False, overwrite_a=False, check_finite=False)  # dsyrk fills the upper triangle
        x = cho_solve(cf, r[:m], check_finite=False); del cf; out.append(float(yy - r[:m] @ x))
    return out
primes = [p for p in range(2, M + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
fam = {'N\\{1}': list(range(2, M + 1)), 'N\\{2}': [k for k in range(1, M + 1) if k != 2],
       'odd': list(range(1, M + 1, 2)), '{1}+primes': [1] + primes}
for sd in range(3):
    rng = np.random.default_rng(100 + sd)
    fam[f'rand1/2 #{sd}'] = [1] + [k for k in range(2, M + 1) if rng.random() < 0.5]
Ms = [m for m in (300, 1000, 1500, 2000, 3000, 4000) if m <= M]
print('family  eta_Mc^2 at Mc = ' + ', '.join(map(str, Ms)), flush=True)
for nm, S in fam.items():
    print(f'{nm:14s} ' + ' '.join(f'{v:.7e}' for v in floors(S, Ms)) + f'   (rss {rss():.0f} MB)', flush=True)
for m in (997, 1001, 1999, 2002, 2999):
    if m <= M and mu[m] != 0:
        v = floors([k for k in range(1, m + 1) if k != m], [m])[0]; print(f'N\\{{{m}}} mu={mu[m]:+d} eta2={v:.6e} m^2*eta2={m*m*v:.5f}', flush=True)
print(f'done {time.time()-t0:.1f}s peak rss {rss():.0f} MB', flush=True)
