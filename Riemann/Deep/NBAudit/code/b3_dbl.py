# nbaudit b3: independent DOUBLE-precision sweep, own code (full cot sums, no symmetry trick, no fixed point): d_N^2 (V1) for all N <= 1280,
# V2 (constraint sum c_k/k = 0) for all n <= 1280; signs of d_N^2 log N - C; comparisons with the s2-copy certified intervals (N <= 640), the shipped
# table (NBRate REPORT), 2611.03924 Table tab:dn, Zenodo 21898513 (d_300^2 interval, cache_gram_300.npy). NOT rigorous (IEEE double).
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 b3_dbl.py 1280
import sys, math, time
import numpy as np, scipy.linalg as sla
from math import gcd
t0 = time.time(); NM = int(sys.argv[1]) if len(sys.argv) > 1 else 1280
LG = math.log(2 * math.pi) - 0.57721566490153286061
import resource
V = {}   # V[(a,k)] = sum_{m=1}^{k-1} {m a/k} cot(pi m/k), a coprime to k, 1 <= a < k
bufI = np.empty(NM * NM, dtype=np.int64); bufF = np.empty(NM * NM)   # preallocated: per-k temporaries made RSS grow to 711 MB (macOS malloc)
for k in range(2, NM + 1):
    m = np.arange(1, k); ct = 1.0 / np.tan(np.pi * m / k)
    a = np.array([x for x in range(1, k) if gcd(x, k) == 1])
    sI = bufI[:len(a) * (k - 1)].reshape(len(a), k - 1); sF = bufF[:len(a) * (k - 1)].reshape(len(a), k - 1)
    np.multiply.outer(a, m, out=sI); np.remainder(sI, k, out=sI); sF[...] = sI
    vals = (sF @ ct) / k
    for x, v in zip(a.tolist(), vals.tolist()): V[(x, k)] = v
    if k % 256 == 0: print(k, 'maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, f'{time.time()-t0:.1f}s', flush=True)
G = np.zeros((NM, NM))
for kk in range(1, NM + 1):
    for hh in range(1, kk + 1):
        if gcd(hh, kk) != 1: continue
        g = LG if hh == kk else (LG * (hh + kk) + (kk - hh) * math.log(hh / kk) - math.pi * (V[(hh, kk)] + (V[(kk % hh, hh)] if hh > 1 else 0.0))) / (2 * hh * kk)
        for d in range(1, NM // kk + 1): G[d * hh - 1, d * kk - 1] = G[d * kk - 1, d * hh - 1] = g / d
del V
ks = np.arange(1, NM + 1); b = (np.log(ks) + 1 - 0.57721566490153286061) / ks
print(f'[{time.time()-t0:.1f}s] Gram built', flush=True)
Lc = np.linalg.cholesky(G); y = sla.solve_triangular(Lc, b, lower=True); del Lc
d1 = np.concatenate([[np.nan], 1 - np.cumsum(y * y)])          # d1[N]
C = 2 + 0.57721566490153286061 - math.log(4 * math.pi)
s = np.sign(d1[2:] * np.log(np.arange(2, NM + 1)) - C)
runs = []; i = 0
while i < len(s):
    j = i
    while j + 1 < len(s) and s[j + 1] == s[i]: j += 1
    runs.append((int(s[i]), i + 2, j + 2)); i = j + 1
print('V1 double sign runs on [2,NM]:', runs, ' sign changes:', len(runs) - 1)
marg = np.abs(d1[2:] * np.log(np.arange(2, NM + 1)) - C)
for lo_, hi_ in ((2, min(640, NM)), (641, NM)):
    if hi_ < lo_: continue
    sl = slice(lo_ - 2, hi_ - 1); k0 = int(np.argmin(marg[sl])) + lo_
    print(f'min |d^2 logN - C| on [{lo_},{hi_}] = {marg[sl].min():.3e} at N = {k0}')
from fractions import Fraction as Fr
mid = {}
for line in open('/tmp/nbaudit/loh_V1.txt'):
    N, a1, b1, c1, e1 = map(int, line.split()); mid[N] = float((Fr(a1, b1) + Fr(c1, e1)) / 2)
dif = np.array([abs(d1[N] - mid[N]) for N in range(1, min(NM, 640) + 1)])
print(f'max |double - certified mid| N<=640: {dif.max():.3e} (at N={int(np.argmax(dif))+1}); N<=100: {dif[:100].max():.3e}')
tab = {695: 0.0067526049812340, 1280: 0.0063611026375655, 640: 0.0070170593529348, 584: 0.0072490880608171}
print('shipped table vs double:', {N: f'{d1[N]-v:+.2e}' for N, v in tab.items() if N <= NM})
dd = np.diff(d1[1:])
print('V1 monotone (d_{N+1}^2 <= d_N^2 + 1e-13) for all N:', bool((dd <= 1e-13).all()), ' max increase', dd.max())
f = d1[2:] * np.log(np.arange(2, NM + 1)); d2f = np.abs(np.diff(f, 2))
print(f'|2nd diff of d^2 logN|: max {d2f.max():.2e} at N={int(np.argmax(d2f))+3}; on (640,NM]: max {d2f[638:].max() if NM > 641 else 0:.2e}; median {np.median(d2f):.2e}')
# V2
kk = np.arange(2, NM + 1)
A = G[1:, 1:] - np.outer(G[1:, 0], 1 / kk) - np.outer(1 / kk, G[0, 1:]) + G[0, 0] / np.outer(kk, kk)
b2 = b[1:] - b[0] / kk
del G
LA = np.linalg.cholesky(A); y2 = sla.solve_triangular(LA, b2, lower=True); del LA
d2 = np.concatenate([[np.nan, np.nan], 1 - np.cumsum(y2 * y2)])  # d2[n]
for n, ref in ((40, '1.274e-2 / s3 0.01273583539133112'), (160, '9.357e-3 / s3 0.009356604331758195'), (640, '7.023e-3 / s3 0.007023264380134277'), (1280, '6.366e-3 / s3 0.006366197363433226'), (320, 's3 0.007961688783869451'), (80, 's3 0.01099879595104947')):
    if n <= NM: print(f'V2 n={n}: {d2[n]:.15e}  (2611.03924 / NBRate s3: {ref})')
if NM >= 300:
    print(f'V2 n=300: {d2[300]:.12e} vs Zenodo 21898513 interval [8.00293986e-3, 8.00294005e-3]; V1 N=300: {d1[300]:.12e}')
    Z = np.load('/Users/adriansager/playground/riemann/.scratch/nbaudit/z21898513/zenodo_package/data/cache_gram_300.npy')
    print('Zenodo cache_gram_300 shape', Z.shape, 'Z[0,0]', Z[0, 0], 'Z[0,1]', Z[0, 1])
    for off in (0, 1):
        n = min(Z.shape[0] - off, 299)
        try: print(f'  offset {off}: max |A - Z| = {np.abs(A[:n, :n] - Z[off:off + n, off:off + n]).max():.3e}')
        except Exception as e: print('  offset', off, e)
np.save('/tmp/nbaudit/d1_b3.npy', d1); np.save('/tmp/nbaudit/d2_b3.npy', d2)
print(f'[{time.time()-t0:.1f}s] done')
