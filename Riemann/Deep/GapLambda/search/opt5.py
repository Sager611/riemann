"""Untrusted optimiser on the well-conditioned FastL tensors (P in shifted-Legendre basis).
usage: opt5.py r dmax thn thd nstarts"""
import sys, time, json, os, resource
sys.path.insert(0, '.')
import numpy as np
from flint import fmpq
import fast, fast2
from opt3 import bfgs


def get_T(r, dmax, thn, thd):
    fn = f'TL_r{r}_d{dmax}_th{thn}_{thd}.npz'
    if os.path.exists(fn):
        d = np.load(fn); return {k: d[k] for k in d.files}
    F = fast2.FastL(fmpq(thn, thd), r, dmax); T = F.tensors(); np.savez(fn, **T); return T


GUARD = 1e-9


def kap2L(T, pl, u, v):
    """kappa^2; returns -1 where float cancellation could matter (c or c1 below GUARD x absolute-value bound)."""
    if v < 0: return -1.0
    ap = np.abs(pl)
    c = T['c'][0] @ pl @ pl @ pl @ pl
    c1v = T['c1'] @ pl @ pl @ pl @ pl
    c1 = np.polyval(c1v[::-1], u)
    c2 = np.polyval((T['c2'] @ pl @ pl @ pl @ pl)[::-1], u)
    ca = np.abs(T['c'][0]) @ ap @ ap @ ap @ ap
    c1a = np.polyval((np.abs(T['c1']) @ ap @ ap @ ap @ ap)[::-1], abs(u))
    if c1 <= GUARD * c1a or c <= GUARD * ca: return -1.0
    return fast.kappa2(c, c1, c2, v)


if __name__ == '__main__':
    r, dmax, thn, thd, ns = map(int, sys.argv[1:6])
    t0 = time.time(); TT = get_T(r, dmax, thn, thd)
    print('tensors %.1fs  peak RSS %.0f MB' % (time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20), flush=True)
    rng = np.random.default_rng(1000 * r + dmax)
    res = {}; prev = None
    for d in range(0, dmax + 1):
        T = {k: A[(slice(None),) + (slice(0, d + 1),) * 4] for k, A in TT.items()}
        f = lambda x: -kap2L(T, x[:d + 1], x[d + 1], x[d + 2]) + (x[:d + 1] @ x[:d + 1] - 1.0) ** 2
        starts = [np.concatenate([[1.0], np.zeros(d), [0.6, 0.4]])]
        if prev is not None:
            starts.append(np.concatenate([prev[:d], [0.0], prev[d:]]))
        for _ in range(ns):
            p = rng.normal(0, 1, d + 1) / np.sqrt(1 + np.arange(d + 1)); p /= np.linalg.norm(p)
            starts.append(np.concatenate([p, [rng.uniform(0.4, 0.9), rng.uniform(0.2, 0.7)]]))
        best = None; vals = []
        for s in starts:
            if f(s) >= 0: continue
            x, fx = bfgs(f, s); x, fx = bfgs(f, x)
            lam = np.sqrt(max(kap2L(T, x[:d + 1], x[d + 1], x[d + 2]), 0)); vals.append(lam)
            if best is None or lam > best[1]: best = (x, lam)
        prev = best[0]; x = best[0]
        pl = x[:d + 1] / x[0]
        res[d] = {'P_leg': list(pl), 'u': x[d + 1], 'v': x[d + 2], 'Lambda': best[1]}
        vals = np.sort(vals)[::-1]
        print('r=%d d=%d Lambda %.6f u=%.5f v=%.5f pL=%s | top %s (%.0fs)' % (r, d, best[1], x[d + 1], x[d + 2],
              np.array2string(pl, precision=4, max_line_width=250), np.array2string(vals[:4], precision=5), time.time() - t0), flush=True)
    json.dump(res, open(f'opt5_r{r}_d{dmax}_th{thn}_{thd}.json', 'w'), indent=1)
    print('peak RSS %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
