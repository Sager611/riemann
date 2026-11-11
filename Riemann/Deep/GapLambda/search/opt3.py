"""Untrusted optimiser: maximise kappa^2(P,u,v) (keyineq3, TeX l.674-676) over P (shifted-Legendre basis, p0=1), u, v.
usage: opt3.py r dmax thn thd [nstarts]"""
import sys, time, json, os
sys.path.insert(0, '.')
import numpy as np
from numpy.polynomial import Legendre, Polynomial
from flint import fmpq
import fast


def tensors(r, dmax, thn, thd):
    fn = f'T_r{r}_d{dmax}_th{thn}_{thd}.npz'
    if os.path.exists(fn):
        d = np.load(fn); return {k: d[k] for k in d.files}
    F = fast.Fast(fmpq(thn, thd), r, dmax)
    T = F.tensors(); np.savez(fn, **T); return T


def legmat(D):
    B = np.zeros((D, D))
    for n in range(D):
        c = Legendre.basis(n, domain=[0, 1]).convert(kind=Polynomial).coef
        B[:len(c), n] = c
    return B


def cvals(T, P):
    return {k: A @ P @ P @ P @ P for k, A in T.items()}


def kap2(T, P, u, v):
    if v < 0: return -1.0
    cs = cvals(T, P)
    c = cs['c'][0]; c1 = np.polyval(cs['c1'][::-1], u); c2 = np.polyval(cs['c2'][::-1], u)
    if c1 <= 0 or c <= 0: return -1.0
    return fast.kappa2(c, c1, c2, v)


def bfgs(f, x0, iters=400, gtol=1e-11):
    n = len(x0)
    def grad(x):
        g = np.zeros(n)
        for i in range(n):
            h = 1e-6 * max(1.0, abs(x[i])); e = np.zeros(n); e[i] = h
            g[i] = (f(x + e) - f(x - e)) / (2 * h)
        return g
    x = np.array(x0, float); fx = f(x); g = grad(x); H = np.eye(n)
    for it in range(iters):
        d = -H @ g
        if g @ d >= 0: H = np.eye(n); d = -g
        t = 1.0
        while True:
            xn = x + t * d; fn = f(xn)
            if fn <= fx + 1e-4 * t * (g @ d) or t < 1e-10: break
            t *= 0.5
        if t < 1e-10 or fx - fn < 1e-15 * abs(fx): 
            if fn < fx: x, fx = xn, fn
            break
        gn = grad(xn); s = xn - x; y = gn - g
        if y @ s > 1e-16:
            rho = 1.0 / (y @ s); I = np.eye(n)
            H = (I - rho * np.outer(s, y)) @ H @ (I - rho * np.outer(y, s)) + rho * np.outer(s, s)
        x, fx, g = xn, fn, gn
        if np.linalg.norm(g) < gtol: break
    return x, fx


if __name__ == '__main__':
    r, dmax, thn, thd = map(int, sys.argv[1:5]); ns = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    t0 = time.time(); T = tensors(r, dmax, thn, thd); D = dmax + 1
    print('tensors %.1fs' % (time.time() - t0), flush=True)
    B = legmat(D); rng = np.random.default_rng(r * 100 + dmax)
    def unpack(x, d):
        p = np.zeros(D); p[:d + 1] = x[:d + 1]
        return B @ p, x[d + 1], x[d + 2]
    res = {}; best_prev = None
    L1 = kap2(T, np.eye(D)[0], 0.6, 0.4)
    print('r=%d P=1 u=.6 v=.4: Lambda %.6f' % (r, np.sqrt(L1)), flush=True)
    for d in range(0, dmax + 1):
        f = lambda x: -kap2(T, *unpack(x, d)) + (x[:d + 1] @ x[:d + 1] - 1.0) ** 2
        starts = []
        if best_prev is not None: starts.append(np.concatenate([best_prev[:d], [0.0], best_prev[d:]]))
        starts.append(np.concatenate([[1.0], np.zeros(d), [0.6, 0.4]]))
        for _ in range(ns):
            p = rng.normal(0, 1, d + 1); p /= np.linalg.norm(p)
            starts.append(np.concatenate([p, [rng.uniform(0.3, 1.0), rng.uniform(0.2, 0.8)]]))
        best = None; vals = []
        for s in starts:
            if f(s) >= 0: continue
            x, fx = bfgs(f, s)
            x, fx = bfgs(f, x)
            vals.append(-fx)
            if best is None or fx < best[1]: best = (x, fx)
        vals = np.sort(np.sqrt(np.maximum(vals, 0)))[::-1]
        print('   top local optima:', np.array2string(vals[:6], precision=5), flush=True)
        best_prev = best[0]
        P, u, v = unpack(best[0], d); P = P / P[0]; Lam = float(np.sqrt(kap2(T, P, u, v)))
        res[d] = {'P_mono': list(P[:d + 1]), 'u': u, 'v': v, 'Lambda': Lam}
        print('r=%d d=%d Lambda %.6f u=%.5f v=%.5f P=%s (%.0fs)' % (r, d, Lam, u, v,
              np.array2string(P[:d + 1], precision=4, max_line_width=200), time.time() - t0), flush=True)
    json.dump(res, open(f'opt3_r{r}_d{dmax}_th{thn}_{thd}.json', 'w'), indent=1)
