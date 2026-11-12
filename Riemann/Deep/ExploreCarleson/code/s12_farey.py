# explore42 s12: lean Farey relaxation of the continuous Nyman-Beurling distance D(1/N)^2 (zeta), Bd large.
# Dilations theta = a/b in [1/N, 1], b <= Bd, plus 1/k (k <= N). Gram <A(al/t),A(be/t)> = a c g(b c, a d) (explore38 s10 / explore13 farey.py,
# Vasyunin's formula). Unique reduced keys computed once (grouped by modulus, cot table per modulus); row-block assembly; Cholesky solve.
# usage: uv run --no-project --with numpy --with scipy python3 s12_farey.py N Bd1,Bd2,...
import sys, math, time, resource
import numpy as np
from math import gcd
from scipy.linalg import cho_factor, cho_solve, eigh
EG = 0.57721566490153286061
N = int(sys.argv[1]); BDS = [int(x) for x in sys.argv[2].split(',')]
def rss(tag):
    if DBG: print('  rss', tag, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6, flush=True)
DBG = len(sys.argv) > 3
def farey(N, Bd):
    s = {(a // gcd(a, b), b // gcd(a, b)) for b in range(1, Bd + 1) for a in range(1, b + 1) if a * N >= b} | {(1, k) for k in range(1, N + 1)}
    return sorted(s, key=lambda t: t[0] / t[1])
def vsums(keys_h, keys_k):
    # V(h,k) = sum_{m=1}^{k-1} {m h/k} cot(pi m/k) for arrays (h,k) with gcd 1; grouped by k
    out = np.zeros(len(keys_h))
    order = np.argsort(keys_k, kind='stable'); kk = keys_k[order]; hh = keys_h[order]
    starts = np.flatnonzero(np.r_[True, kk[1:] != kk[:-1]]); ends = np.r_[starts[1:], len(kk)]
    for s0, e0 in zip(starts, ends):
        k = int(kk[s0])
        if k == 1: continue
        m = np.arange(1, k, dtype=np.int64); ct = 1.0 / np.tan(np.pi * m / k)
        hs = hh[s0:e0]
        for c0 in range(0, len(hs), max(1, 250000 // k)):
            hc = hs[c0:c0 + max(1, 250000 // k)]
            r = (hc[:, None] * m[None, :]) % k
            out[order[s0 + c0:s0 + c0 + len(hc)]] = (r * ct[None, :]).sum(axis=1) / k
    return out
def run(N, Bd):
    t0 = time.time()
    fr = farey(N, Bd); n = len(fr)
    A = np.array([a for a, b in fr], dtype=np.int64); B = np.array([b for a, b in fr], dtype=np.int64)
    # pass 1: unique reduced keys over i <= j
    allkeys = []
    for i0 in range(0, n, 256):
        i1 = min(n, i0 + 256)
        h = B[i0:i1, None] * A[None, :]; k = A[i0:i1, None] * B[None, :]
        g = np.gcd(h, k); h //= g; k //= g
        lo = np.minimum(h, k); hi = np.maximum(h, k)
        allkeys.append(np.unique(lo * 1000003 + hi))
        allkeys = [np.unique(np.concatenate(allkeys))]
    keys = allkeys[0]; lo = keys // 1000003; hi = keys % 1000003
    rss('pass1')
    vlh = vsums(lo, hi); vhl = vsums(hi, lo)   # V(lo,hi), V(hi,lo)
    rss('vsums')
    lof = lo.astype(float); hif = hi.astype(float)
    gzv = (math.log(2 * math.pi) - EG) * (lof + hif) / (2 * lof * hif) + (hif - lof) / (2 * lof * hif) * np.log(lof / hif) - math.pi / (2 * lof * hif) * (vlh + vhl)
    G = np.empty((n, n), order='F')
    for i0 in range(0, n, 256):
        i1 = min(n, i0 + 256)
        h = B[i0:i1, None] * A[None, :]; k = A[i0:i1, None] * B[None, :]
        g = np.gcd(h, k); h //= g; k //= g
        kl = np.minimum(h, k) * 1000003 + np.maximum(h, k)
        val = gzv[np.searchsorted(keys, kl)]
        G[i0:i1, :] = (A[i0:i1, None] * A[None, :]) * val / g
    rss('G')
    th = A / B
    bb = th * (np.log(1 / th) + 1 - EG)
    w = np.linalg.eigvalsh(G) if n <= 1500 else None
    rss('eig')
    cf = cho_factor(G, lower=True, overwrite_a=True, check_finite=False)
    y = cho_solve(cf, bb, check_finite=False)
    d2 = 1 - float(bb @ y)
    L = math.log(N)
    cinv = (w.min() / w.max()) if w is not None else float('nan')
    print(f'N={N} Bd={Bd} n={n} keys={len(keys)} d2={d2:.10f} d2logN={d2*L:.8f} cond^-1={cinv:.2e} {time.time()-t0:.1f}s maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB', flush=True)
for Bd in BDS:
    run(N, Bd)
