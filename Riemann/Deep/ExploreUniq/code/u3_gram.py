# explore37 u3: actual Baez-Duarte V1 distances d_N(S)^2 for dilation sets S (Vasyunin Gram, double, copied from
# Riemann/Deep/NBRate/code/s1_double.py) against the unconditional time-domain floors eta_M(S)^2 of u1_exact.py;
# minimality constants delta_{m,N}^2 = 1/(G_N^{-1})_{mm}; optimal coefficients c*_k(N) vs mu(k).
# usage: uv run --no-project --with numpy python3 u3_gram.py NMAX MTIME
import sys, math, time, resource, json
import numpy as np
from math import gcd
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 640
MT = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
EG = 0.57721566490153286061
C = 2 + EG - math.log(4 * math.pi)
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
V = np.zeros((NMAX + 1, NMAX + 1))
for k in range(2, NMAX + 1):
    h = (k - 1) // 2
    if h == 0: continue
    m = np.arange(1, h + 1)
    a = np.array([x for x in range(1, k) if gcd(x, k) == 1])
    R = 2 * (np.outer(a, m) % k) - k
    V[a, k] = (R @ (1.0 / np.tan(np.pi * m / k))) / k
LG = math.log(2 * math.pi) - EG
def g(h, k):
    vhk = V[h, k] if k > 1 else 0.0
    vkh = V[k % h, h] if h > 1 else 0.0
    return (LG * (h + k) + (k - h) * math.log(h / k) - math.pi * (vhk + vkh)) / (2 * h * k)
G = np.zeros((NMAX, NMAX))
for hh in range(1, NMAX + 1):
    for kk in range(hh, NMAX + 1):
        d = gcd(hh, kk)
        G[hh - 1, kk - 1] = G[kk - 1, hh - 1] = g(hh // d, kk // d) / d
del V
ks = np.arange(1, NMAX + 1)
b = (np.log(ks) + 1 - EG) / ks
print(f'Gram {NMAX} built {time.time()-t0:.1f}s rss {rss():.0f} MB; G11={G[0,0]:.12f}', flush=True)
mu = np.ones(NMAX + 1, dtype=np.int64); mu[0] = 0
for p in range(2, NMAX + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu[p::p] *= -1; mu[p * p::p * p] = 0
def dS(S, N):
    idx = [k - 1 for k in S if k <= N]
    Gs = G[np.ix_(idx, idx)]; bs = b[idx]
    L = np.linalg.cholesky(Gs); y = np.linalg.solve(L, bs)
    return 1 - y @ y
Ns = [n for n in (20, 40, 80, 160, 320, 640, 1280) if n <= NMAX]
# full d_N^2 and optimal coefficients
print('N  d_N^2  d_N^2 logN   -c*_k(N)/mu(k) for k=1,2,3,5,6,7,10,11 ; max|c*_k| over non-squarefree k<=10', flush=True)
cstar = {}
for N in Ns:
    c = np.linalg.solve(G[:N, :N], b[:N]); cstar[N] = c
    d2 = 1 - b[:N] @ c
    rat = [-c[k - 1] / mu[k] for k in (1, 2, 3, 5, 6, 7, 10, 11)]  # chi = -sum mu(k) rho_{1/k} (Moebius identity)
    nsq = max(abs(c[k - 1]) for k in (4, 8, 9) if k <= N)
    print(f'{N:5d} {d2:.10f} {d2*math.log(N):.6f}  ' + ' '.join(f'{r:+.4f}' for r in rat) + f'  {nsq:.4f}', flush=True)
# minimality constants and S = N\{m}
print('m  mu  delta2_{m,N} (N=' + ','.join(map(str, Ns)) + ')  | d_N(N\\{m})^2 at same N | m^2*delta2 at N=max', flush=True)
for m in (1, 2, 3, 4, 5, 6, 7, 10, 30, 101):
    dl = []; dd = []
    for N in Ns:
        if N < m: dl.append(float('nan')); dd.append(float('nan')); continue
        Gi = np.linalg.inv(G[:N, :N]); dl.append(1 / Gi[m - 1, m - 1])
        c = cstar[N]; d2 = 1 - b[:N] @ c
        dd.append(d2 + c[m - 1] ** 2 / Gi[m - 1, m - 1])
    print(f'{m:4d} {mu[m]:+d} ' + ' '.join(f'{x:.5e}' for x in dl) + ' | ' + ' '.join(f'{x:.5e}' for x in dd) + f' | {m*m*dl[-1]:.4f}', flush=True)
# families: same definitions as u1 but random sets from per-seed generators; floors recomputed here (time domain, MT)
primes = [p for p in range(2, max(NMAX, MT) + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
top = max(NMAX, MT)
fam = {'N\\{1}': list(range(2, top + 1)), 'N\\{2}': [k for k in range(1, top + 1) if k != 2],
       'odd': list(range(1, top + 1, 2)), '{1}+primes': [1] + primes,
       'squarefree': [k for k in range(1, top + 1) if mu[min(k, NMAX)] != 0 or k > NMAX],
       'N\\{4}': [k for k in range(1, top + 1) if k != 4]}
for sd in range(3):
    rng = np.random.default_rng(100 + sd)
    fam[f'rand1/2 #{sd}'] = [1] + [k for k in range(2, top + 1) if rng.random() < 0.5]
# time-domain floors (same construction as u1_exact.py)
import mpmath as mp
mp.mp.dps = 40
M = MT
R00 = np.zeros(M + 1); R01 = np.zeros(M + 1); R11 = np.zeros(M + 1)
for j in range(1, M + 1):
    a_ = mp.mpf(1) / (j + 1); b_ = mp.mpf(1) / j; Lj = mp.log(j)
    F1 = lambda x: x * (mp.log(x) + Lj) - x
    F2 = lambda x: (lambda lx: x * lx**2 - 2 * x * lx + 2 * x)(mp.log(x) + Lj)
    m0 = b_ - a_; m1 = F1(b_) - F1(a_); m2 = F2(b_) - F2(a_)
    r00 = mp.sqrt(m0); r01 = m1 / r00; r11 = mp.sqrt(m2 - r01**2)
    R00[j], R01[j], R11[j] = float(r00), float(r01), float(r11)
js = np.arange(1, M + 1); lj = np.log(js.astype(float))
def block(alpha, beta, jj): return np.stack([R00[jj] * alpha + R01[jj] * beta, R11[jj] * beta], axis=1).reshape(-1)
evec = block(1 - lj, np.ones(M), js)
VT = np.zeros((2 * M, M)); lgam = np.array([math.lgamma(n + 1) for n in range(M + 2)])
for k in range(1, M + 1):
    jj = np.arange(k, M + 1); n = jj // k
    VT[2 * (k - 1):, k - 1] = block(n * (1 + math.log(k) - np.log(jj.astype(float))) + lgam[n], n.astype(float), jj)
def eta2(S, Mc):
    rows = 2 * Mc; cols = [k - 1 for k in S if k <= Mc]
    A = VT[:rows, cols]; y = evec[:rows]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None); r = y - A @ coef
    return float(r @ r)
print(f'time-domain V {VT.shape} built {time.time()-t0:.1f}s rss {rss():.0f} MB', flush=True)
print('family | floor eta_Mc^2 (Mc=' + ','.join(str(x) for x in (10, 100, M)) + ') | actual d_N(S)^2 (N=' + ','.join(map(str, Ns)) + ')', flush=True)
res = {}
for nm, S in fam.items():
    fl = [eta2(S, Mc) for Mc in (10, 100, M)]
    act = [dS(S, N) for N in Ns]
    res[nm] = {'floor': fl, 'actual': act}
    print(f'{nm:14s} | ' + ' '.join(f'{x:.5e}' for x in fl) + ' | ' + ' '.join(f'{x:.5e}' for x in act), flush=True)
json.dump({'Ns': Ns, 'fam': res}, open('/Users/adriansager/playground/riemann/.scratch/explore37/out/u3_%d_%d.json' % (NMAX, MT), 'w'))
print(f'done {time.time()-t0:.1f}s peak rss {rss():.0f} MB', flush=True)
