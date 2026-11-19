# nbaudit b4: BCF ramp check (COROLLARIES l.51 dependency), double precision, own Gram code (as b3):
# I_N = Q_N(c), c_k = -mu(k)(1 - log k/log N); compare with a2's certified bcf_lo/bcf_hi (N <= 640, verbatim s2 copy) and the window means of
# log^2 N * I_N - C log N quoted by NBRate/COROLLARIES (3.5480 [200,400], 3.5496 [400,800], 3.5478 [800,1280]).
# usage: uv run --no-project --with numpy python3 b4_ramp.py 1280
import sys, math, time, resource
import numpy as np
from math import gcd
t0 = time.time(); NM = int(sys.argv[1]) if len(sys.argv) > 1 else 1280
EG = 0.57721566490153286061; LG = math.log(2 * math.pi) - EG
V = {}; bufI = np.empty(NM * NM, dtype=np.int64); bufF = np.empty(NM * NM)
for k in range(2, NM + 1):
    m = np.arange(1, k); ct = 1.0 / np.tan(np.pi * m / k)
    a = np.array([x for x in range(1, k) if gcd(x, k) == 1])
    sI = bufI[:len(a) * (k - 1)].reshape(len(a), k - 1); sF = bufF[:len(a) * (k - 1)].reshape(len(a), k - 1)
    np.multiply.outer(a, m, out=sI); np.remainder(sI, k, out=sI); sF[...] = sI
    for x, v in zip(a.tolist(), ((sF @ ct) / k).tolist()): V[(x, k)] = v
del bufI, bufF
G = np.zeros((NM, NM))
for kk in range(1, NM + 1):
    for hh in range(1, kk + 1):
        if gcd(hh, kk) != 1: continue
        g = LG if hh == kk else (LG * (hh + kk) + (kk - hh) * math.log(hh / kk) - math.pi * (V[(hh, kk)] + (V[(kk % hh, hh)] if hh > 1 else 0.0))) / (2 * hh * kk)
        for d in range(1, NM // kk + 1): G[d * hh - 1, d * kk - 1] = G[d * kk - 1, d * hh - 1] = g / d
del V
ks = np.arange(1, NM + 1); b = (np.log(ks) + 1 - EG) / ks
mu = np.ones(NM + 1, dtype=np.int64); mu[0] = 0
for p in range(2, NM + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)): mu[p::p] *= -1; mu[p * p::p * p] = 0
C = 2 + EG - math.log(4 * math.pi)
I = np.full(NM + 1, np.nan)
for N in range(2, NM + 1):
    c = -mu[1:N + 1] * (1 - np.log(ks[:N]) / math.log(N))
    I[N] = 1 - 2 * b[:N] @ c + c @ (G[:N, :N] @ c)
print(f'[{time.time()-t0:.1f}s] ramp done; maxrss MB {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f}', flush=True)
z = np.load('/tmp/nbaudit/rig_V1.npz'); lo, hi = z['bcf_lo'], z['bcf_hi']; n0 = min(640, NM)
ins = sum(1 for N in range(2, n0 + 1) if lo[N] - 1e-13 <= I[N] <= hi[N] + 1e-13)
print(f'double ramp within 1e-13 of the certified a2 ramp bounds for {ins} of {n0 - 1} N <= {n0}; max |double - certified mid| = {max(abs(I[N] - (lo[N] + hi[N]) / 2) for N in range(2, n0 + 1)):.2e}')
F = np.array([math.log(N) ** 2 * I[N] - C * math.log(N) if N >= 2 else np.nan for N in range(NM + 1)])
for a_, b_ in ((200, 400), (400, 800), (800, NM)):
    if b_ <= NM: print(f'mean of log^2N I_N - C logN on [{a_},{b_}] = {np.mean(F[a_:b_ + 1]):.4f} (sd {np.std(F[a_:b_ + 1]):.4f})')
if NM >= 1280:
    Lc = np.linalg.cholesky(G); y = np.linalg.solve(Lc, b); d1280 = 1 - y @ y
    print(f'I_1280 = {I[1280]:.4f}, d_1280^2 = {d1280:.10f}, ratio {I[1280] / d1280:.2f}')
