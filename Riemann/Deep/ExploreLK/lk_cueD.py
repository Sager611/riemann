#!/usr/bin/env python3
"""explore58: add Gonek's shifted mean values D(a) = <sum_gamma Z(gamma+a)^2> / <Z^2> = 1 - (sin pi a / pi a)^2 (RH; normalised units)
to the cone test. Rows: m_0..m_J (pi^{2j}/(2j+1)), d = pi^2/3, D(a) on a grid. Columns: (hard-core) CUE char polys, period N, scaled by s.
"""
import sys, math, json
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, 'Riemann/Deep/ExploreLK')
from lk_cue import cue_angles, hardcore
PI = math.pi

def fvals(xs, P, pts):
    lf = np.sum(np.log(np.abs(np.sin(PI * (pts[:, None] - xs[None, :]) / P)) + 1e-300), axis=1)
    sg = np.prod(np.sign(np.sin(PI * (pts[:, None] - xs[None, :]) / P)), axis=1)
    return sg, lf

def column(xs, P, s, J, alphas):
    n2 = len(xs); n = n2 // 2; M = 4 * n2 + 8
    t = np.arange(M) * (P / M)
    sg, lf = fvals(xs, P, t); sh = np.max(lf); f = sg * np.exp(lf - sh)
    c = np.fft.fft(f) / M; p = np.abs(c) ** 2
    k = np.rint(np.fft.fftfreq(M, d=1.0 / M)).astype(int)
    w = np.zeros(n + 1)
    for kk, pp in zip(k, p):
        if abs(kk) <= n: w[abs(kk)] += pp
    m0 = w.sum(); om = s * 2 * PI * np.arange(n + 1) / P
    mj = [np.sum(w * om ** (2 * j)) / m0 for j in range(J + 1)]
    S = np.sin(PI * (xs[:, None] - xs[None, :]) / P); np.fill_diagonal(S, 1.0)
    ls = np.sum(np.log(np.abs(S)), axis=1) - sh
    d = np.sum((PI / P) ** 2 * np.exp(2 * ls)) / P / m0 * s ** 3
    Ds = []
    for a in alphas:  # g(x) = f(s x): zeros x_k/s; D_g(a) = (s/P) sum_k f(x_k + s a)^2 / m0
        sgv, lv = fvals(xs, P, xs + s * a)
        Ds.append(s / P * np.sum(np.exp(2 * (lv - sh))) / m0)
    return np.array(mj + [d] + Ds)

def target(J, alphas):
    return np.array([PI ** (2 * j) / (2 * j + 1) for j in range(J + 1)] + [PI ** 2 / 3] + [1 - (math.sin(PI * a) / (PI * a)) ** 2 for a in alphas])

def lp(V, t):
    A = V.T; R = 1.0 / t; A = A * R[:, None]; b = t * R
    nr, nc = A.shape
    res = linprog(np.concatenate([np.zeros(nc), np.ones(2 * nr)]), A_eq=np.hstack([A, np.eye(nr), -np.eye(nr)]), b_eq=b,
                  bounds=[(0, None)] * (nc + 2 * nr), method='highs')
    r = res.x[nc:nc + nr] - res.x[nc + nr:]
    return res.fun / nr, r

if __name__ == '__main__':
    N = int(sys.argv[1]); S = int(sys.argv[2]); deltas = [float(a) for a in sys.argv[3].split(',')]; J = int(sys.argv[4])
    alphas = [float(a) for a in sys.argv[5].split(',')]
    rng = np.random.default_rng(97)
    samples = [cue_angles(N, rng) * N / (2 * PI) for _ in range(S)]
    t = target(J, alphas); out = {}
    for delta in deltas:
        cols = []
        for xs in samples:
            ys = hardcore(xs, float(N), delta + 1e-6)[0] if delta > 0 else xs
            for s in [1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5]:
                cols.append(column(ys, float(N), s, J, alphas))
        cols.append(np.concatenate([[1.0], np.zeros(J), [0.0], np.zeros(len(alphas))]) * 0 + column(np.arange(N) + 0.5, float(N), 1.0, J, alphas))  # sine
        val, r = lp(np.array(cols), t)
        print(f'N={N} delta={delta:.2f} J={J} alphas={alphas}: mean|rel resid|={val:.3e}; rel resid per row={np.round(r,4).tolist()}', flush=True)
        out[str(delta)] = dict(resid=val, rows=r.tolist())
    json.dump(out, open(f'Riemann/Deep/ExploreLK/logs/cueD_N{N}_J{J}.json', 'w'), indent=1)
