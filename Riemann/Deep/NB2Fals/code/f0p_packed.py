# nb2fals f0p: as f0_perN.py but packed storage (upper, column-major, 4N^2 bytes) so that N ~ 5500 fits under the 290 MB cap.
# A = U^T U (dpptrf); the leading N-block of U factors G_N (packed prefix). Per N: v = U_N c_N (dtpmv on the prefix), S_N = L^2||v - y_N||^2.
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 f0p_packed.py NMAX -> .scratch/nb2fals/perN{NMAX}.npz
import sys, math, time, resource
import numpy as np, mpmath as mp
from scipy.linalg import lapack, blas
t0 = time.time()
NMAX = int(sys.argv[1])
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
OUT = '/Users/adriansager/playground/riemann/.scratch/nb2fals/'
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
off = np.zeros(NMAX + 2, dtype=np.int64)
for k in range(1, NMAX + 1): off[k + 1] = off[k] + k // 2 + 1
Vf = np.zeros(int(off[NMAX + 1])); BI = np.empty(256 * (NMAX // 2 + 1), dtype=np.int64); BF = np.empty(256 * (NMAX // 2 + 1))
for k in range(3, NMAX + 1):
    h = (k - 1) // 2; m = np.arange(1, h + 1); ct = 1.0 / np.tan(np.pi * m / k)
    a = m[np.gcd(m, k) == 1]
    for i in range(0, len(a), 256):
        aa = a[i:i + 256]; P = BI[:len(aa) * h].reshape(len(aa), h); np.multiply.outer(aa, m, out=P); np.remainder(P, k, out=P)
        Pf = BF[:len(aa) * h].reshape(len(aa), h); np.multiply(P, 2, out=Pf, casting='unsafe'); Pf -= k
        Vf[off[k] + aa] = (Pf @ ct) / k
t1 = time.time(); print(f'  V built {t1-t0:.0f}s rss {rss():.0f} MB', flush=True)
def Vlook(a, k):
    s = a <= k // 2; idx = off[k] + np.where(s, a, k - a)
    return np.where(s, Vf[idx], -Vf[idx])
LG = math.log(2 * math.pi) - EG
AP = np.zeros(NMAX * (NMAX + 1) // 2)
for kk in range(1, NMAX + 1):            # column kk: rows hh = 1..kk
    hh = np.arange(1, kk + 1); d = np.gcd(hh, kk); h = hh // d; k = kk // d
    vhk = np.zeros(kk); m1 = k > 1; vhk[m1] = Vlook(h[m1], k[m1])
    vkh = np.zeros(kk); m2 = h > 1; vkh[m2] = Vlook(k[m2] % h[m2], h[m2])
    g = (LG * (h + k) + (k - h) * np.log(h / k) - np.pi * (vhk + vkh)) / (2 * h * k)
    j = kk - 1; AP[j * (j + 1) // 2: j * (j + 1) // 2 + kk] = g / d
del Vf, BI, BF
G11, G12, G23 = AP[0], AP[1], AP[4]
t2 = time.time(); print(f'  G built {t2-t1:.0f}s rss {rss():.0f} MB', flush=True)
ks = np.arange(1, NMAX + 1); b = (np.log(ks) + 1 - EG) / ks
mu = np.ones(NMAX + 1, dtype=np.int64); mu[0] = 0
for p in range(2, NMAX + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)): mu[p::p] *= -1; mu[p * p::p * p] = 0
mu = mu[1:].astype(float)
U, info = lapack.dpptrf(NMAX, AP, lower=0, overwrite_ap=1)
assert info == 0 and np.shares_memory(U, AP), (info, np.shares_memory(U, AP))
t3 = time.time(); print(f'  dpptrf {t3-t2:.0f}s rss {rss():.0f} MB', flush=True)
mp.mp.dps = 20
a1 = float(mp.zeta(0, derivative=2) / mp.zeta(0) - (mp.zeta(0, derivative=1) / mp.zeta(0)) ** 2)
z = [0.0] + [float(mp.zeta(2 * k + 1)) for k in range(1, 200)]
def ell(n):
    s0 = -(math.log(2 * math.pi) * (1 + math.log(n) - EG) + a1) / n
    if n == 1: tail = sum((z[k] - 1) / (2 * k * (2 * k + 1)) for k in range(1, 199)) + (1 - math.log(2))
    else:
        x = 1.0 / n; tail = sum(z[k] * x ** (2 * k + 1) / (2 * k * (2 * k + 1)) for k in range(1, 60))
    return s0 - tail
lth = np.array([ell(int(n)) for n in ks])
y = blas.dtpsv(NMAX, U, b.copy(), lower=0, trans=1); w = blas.dtpsv(NMAX, U, lth.copy(), lower=0, trans=1)
d2 = 1 - np.cumsum(y * y); Sig = np.cumsum(w * w)
Ns = np.arange(2, NMAX + 1); out = np.zeros((len(Ns), 8))
for i, N in enumerate(Ns):
    L = math.log(N); cr = -mu[:N] * (1 - np.log(ks[:N]) / L)
    v = blas.dtpmv(N, U, cr.copy(), lower=0, trans=0); Q = float(v @ v); I = 1 - 2 * float(cr @ b[:N]) + Q
    S = L * L * float(((v - y[:N]) ** 2).sum())
    out[i] = (N, L, L * L * d2[N - 1] - C * L, S - Sig[N - 1], J - Sig[N - 1], L * L * I - C * L - J, d2[N - 1], S)
t4 = time.time(); print(f'  perN {t4-t3:.0f}s rss {rss():.0f} MB', flush=True)
rng = np.random.default_rng(1); x = rng.standard_normal(NMAX); x /= np.linalg.norm(x)
for it in range(60):
    x = blas.dtpsv(NMAX, U, x, lower=0, trans=1); x = blas.dtpsv(NMAX, U, x, lower=0, trans=0); lam_inv = np.linalg.norm(x); x /= lam_inv
x = rng.standard_normal(NMAX); x /= np.linalg.norm(x)
for it in range(60):
    x = blas.dtpmv(NMAX, U, x, lower=0, trans=0); x = blas.dtpmv(NMAX, U, x, lower=0, trans=1); lam_max = np.linalg.norm(x); x /= lam_max
chk = np.abs(out[:, 2] - (out[:, 4] + out[:, 5] - out[:, 3])).max()
np.savez(OUT + f'perN{NMAX}.npz', perN=out, cols=np.array(['N', 'L', 'C2', 'kappa', 'JmSig', 'R', 'd2', 'S']), Sig=Sig, d2=d2)
print(f'NMAX={NMAX} (packed): G11={G11:.15f} G12={G12:.15f} G23={G23:.15f}; total {time.time()-t0:.0f}s; maxrss {rss():.0f} MB')
print(f'  lambda_min ~ {1/lam_inv:.4e}, lambda_max ~ {lam_max:.4f}, cond ~ {lam_max*lam_inv:.3e}; identity check max {chk:.2e}')
for NN in (640, 1280, 2560, 4000, 5000, 5500, 6000):
    if NN <= NMAX: print(f'  d2[{NN}] = {d2[NN-1]:.16f}  C2 = {math.log(NN)**2*d2[NN-1]-C*math.log(NN):+.6f}  Sigma = {Sig[NN-1]:.6f}  kappa = {out[NN-2,3]:.6f}')
