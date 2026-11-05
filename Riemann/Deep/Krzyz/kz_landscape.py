#!/usr/bin/env python3
"""Floating-point multistart maximization of |a_n| over configurations with N atoms
(sanity check of the smaller-strata claims sup|a_n| < 2/e; NOT a proof).
usage: kz_landscape.py n N starts seed"""
import sys, math, json
import numpy as np
from scipy.optimize import minimize

def an(x, n, N):
    lam = np.abs(x[:N]); th = np.concatenate([[0.0], x[N:]])
    h = np.zeros(n + 1, dtype=complex)
    for m in range(1, n + 1):
        h[m] = -2 * np.sum(lam * np.exp(1j * m * th))
    f = np.zeros(n + 1, dtype=complex); f[0] = 1
    for k in range(1, n + 1):
        f[k] = sum(j * h[j] * f[k - j] for j in range(1, k + 1)) / k
    return abs(f[n]) * math.exp(-lam.sum())

if __name__ == '__main__':
    n, N, starts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = np.random.default_rng(seed); best = []
    for s in range(starts):
        x0 = np.concatenate([rng.uniform(0.05, 1.5, N), rng.uniform(0, 2 * math.pi, N - 1)])
        r = minimize(lambda x: -an(x, n, N), x0, method='BFGS', options=dict(gtol=1e-10, maxiter=2000))
        best.append((-r.fun, list(np.abs(r.x[:N])), list(np.mod(r.x[N:], 2 * math.pi))))
    best.sort(key=lambda b: -b[0])
    vals = []
    for b in best:
        if all(abs(b[0] - v) > 1e-7 for v in vals): vals.append(b[0])
    two_e = 2 / math.e
    print(json.dumps(dict(n=n, N=N, starts=starts, max=best[0][0], gap_to_2_over_e=two_e - best[0][0], t=sum(best[0][1]),
                          masses=[round(v, 4) for v in best[0][1]], distinct_top=[round(v, 6) for v in vals[:6]])))
