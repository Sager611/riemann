# nb2fals f0: memory-lean double Vasyunin Gram G_N (V1, same formulas as NB2nd t0) for N <= NMAX, in-place Cholesky,
# per-N C2(N) = L^2 d_N^2 - C L, kappa(N) = S_N - Sigma_N, J - Sigma_N, R(N) = L^2 I_N - C L - J  (NB2nd derivation.md eqs (1), (5)).
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 f0_perN.py NMAX  -> .scratch/nb2fals/perN{NMAX}.npz
import sys, math, time, resource
import numpy as np, mpmath as mp
from scipy.linalg import lapack, solve_triangular
t0 = time.time()
NMAX = int(sys.argv[1])
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.5626634565363677
OUT = '/Users/adriansager/playground/riemann/.scratch/nb2fals/'
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
# Vasyunin sums V(a,k) (a coprime to k): stored for a <= k//2 only; V(k-a,k) = -V(a,k) exactly (term-by-term sign flip)
off = np.zeros(NMAX + 2, dtype=np.int64)
for k in range(1, NMAX + 1): off[k + 1] = off[k] + k // 2 + 1
Vf = np.zeros(int(off[NMAX + 1])); BI = np.empty(256 * (NMAX // 2 + 1), dtype=np.int64); BF = np.empty(256 * (NMAX // 2 + 1))
for k in range(3, NMAX + 1):
    h = (k - 1) // 2; m = np.arange(1, h + 1); ct = 1.0 / np.tan(np.pi * m / k)
    a = m[np.gcd(m, k) == 1]
    for i in range(0, len(a), 256):   # fixed work buffers: macOS malloc caches freed large blocks (rss 230 MB at N = 2560 otherwise)
        aa = a[i:i + 256]; P = BI[:len(aa) * h].reshape(len(aa), h); np.multiply.outer(aa, m, out=P); np.remainder(P, k, out=P)
        Pf = BF[:len(aa) * h].reshape(len(aa), h); np.multiply(P, 2, out=Pf, casting='unsafe'); Pf -= k
        Vf[off[k] + aa] = (Pf @ ct) / k
t1 = time.time(); print(f'  rss after V {rss():.0f} MB', flush=True)
def Vlook(a, k):
    s = a <= k // 2; idx = off[k] + np.where(s, a, k - a)
    return np.where(s, Vf[idx], -Vf[idx])
LG = math.log(2 * math.pi) - EG
G = np.zeros((NMAX, NMAX))          # upper triangle filled (row-major); = lower triangle of the F-view G.T
for hh in range(1, NMAX + 1):
    kk = np.arange(hh, NMAX + 1); d = np.gcd(hh, kk); h = hh // d; k = kk // d
    vhk = np.zeros(len(kk)); m1 = k > 1; vhk[m1] = Vlook(h[m1], k[m1])
    vkh = np.zeros(len(kk)); m2 = h > 1; vkh[m2] = Vlook(k[m2] % h[m2], h[m2])
    g = (LG * (h + k) + (k - h) * np.log(h / k) - np.pi * (vhk + vkh)) / (2 * h * k)
    G[hh - 1, hh - 1:] = g / d
del Vf, BI, BF
ks = np.arange(1, NMAX + 1); b = (np.log(ks) + 1 - EG) / ks
mu = np.ones(NMAX + 1, dtype=np.int64); mu[0] = 0
for p in range(2, NMAX + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)): mu[p::p] *= -1; mu[p * p::p * p] = 0
mu = mu[1:].astype(float)
G11, G12, G23 = G[0, 0], G[0, 1], G[1, 2]
# row abs sums of the symmetric G (1-norm) from the upper triangle, chunked
an = np.zeros(NMAX)
for i in range(NMAX): row = np.abs(G[i, i:]); an[i] += row.sum(); an[i + 1:] += row[1:]
anorm = an.max(); t2 = time.time(); print(f'  rss after G {rss():.0f} MB', flush=True)
Fv = G.T                             # F-contiguous view
c, info = lapack.dpotrf(Fv, lower=1, clean=1, overwrite_a=1)
assert info == 0 and np.shares_memory(c, G), (info, np.shares_memory(c, G))
U = G                                # U = L^T (upper, row-major); leading N-block = Cholesky factor of G_N
try: rcond, info2 = lapack.dpocon(c, anorm, uplo='L')
except Exception as ex: rcond = float('nan'); print('dpocon failed', ex)
t3 = time.time(); print(f'  rss after chol {rss():.0f} MB', flush=True)
# head functional l_n (NB2nd eq. (2))
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
y = solve_triangular(c, b, lower=True, check_finite=False); w = solve_triangular(c, lth, lower=True, check_finite=False)
d2 = 1 - np.cumsum(y * y); Sig = np.cumsum(w * w); print(f'  rss after solves {rss():.0f} MB', flush=True)
Ns = np.arange(2, NMAX + 1); out = np.zeros((len(Ns), 8))
for i, N in enumerate(Ns):
    L = math.log(N); cr = -mu[:N] * (1 - np.log(ks[:N]) / L)
    v = U[:N, :N] @ cr; Q = float(v @ v); I = 1 - 2 * float(cr @ b[:N]) + Q
    S = L * L * float(((v - y[:N]) ** 2).sum())
    C2 = L * L * d2[N - 1] - C * L; R = L * L * I - C * L - J
    out[i] = (N, L, C2, S - Sig[N - 1], J - Sig[N - 1], R, d2[N - 1], S)
t4 = time.time(); print(f'  rss after perN {rss():.0f} MB', flush=True)
# lambda_min(G_NMAX) by inverse iteration, lambda_max by power iteration (interlacing: cond(G_N) <= these for N <= NMAX)
rng = np.random.default_rng(1); x = rng.standard_normal(NMAX); x /= np.linalg.norm(x)
for it in range(60):
    x = solve_triangular(c, x, lower=True, check_finite=False); x = solve_triangular(c, x, lower=True, trans='T', check_finite=False)
    lam_inv = np.linalg.norm(x); x /= lam_inv
x = rng.standard_normal(NMAX); x /= np.linalg.norm(x)
for it in range(60):
    x = U.T @ (U @ x); lam_max = np.linalg.norm(x); x /= lam_max
chk = np.abs(out[:, 2] - (out[:, 4] + out[:, 5] - out[:, 3])).max()
np.savez(OUT + f'perN{NMAX}.npz', perN=out, cols=np.array(['N', 'L', 'C2', 'kappa', 'JmSig', 'R', 'd2', 'S']), Sig=Sig, d2=d2)
print(f'NMAX={NMAX}: G11={G11:.15f} G12={G12:.15f} G23={G23:.15f}; times V {t1-t0:.1f}s G {t2-t1:.1f}s chol {t3-t2:.1f}s perN {t4-t3:.1f}s; maxrss {rss():.0f} MB')
print(f'  lambda_min(G_NMAX) ~ {1/lam_inv:.4e}, lambda_max ~ {lam_max:.4f}, cond ~ {lam_max*lam_inv:.3e}; dpocon 1-norm rcond {rcond:.3e}; identity check |C2-(J-Sig+R-kappa)| max {chk:.2e}')
for NN in (640, 695, 1280, 2560, 4000, 5000):
    if NN <= NMAX: print(f'  d2[{NN}] = {d2[NN-1]:.16f}  C2 = {math.log(NN)**2*d2[NN-1]-C*math.log(NN):+.6f}  Sigma = {Sig[NN-1]:.6f}')
for lo, hi in [(100, 250), (250, 450), (450, 640), (800, 1280), (1280, 2560), (2560, 5120)]:
    m = (out[:, 0] >= lo) & (out[:, 0] <= hi)
    if m.sum() > 10 and hi <= NMAX: print(f'  N in [{lo},{hi}]: C2 {out[m,2].mean():+.5f} sd {out[m,2].std():.5f} | kappa {out[m,3].mean():.4f} sd {out[m,3].std():.4f} | J-Sig {out[m,4].mean():.4f} | R {out[m,5].mean():+.4f} sd {out[m,5].std():.4f}')
