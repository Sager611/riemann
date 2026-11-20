#!/usr/bin/env python3
"""ceilaudit7 A2: explore58's cone test (hard-core CUE columns, gaps >= delta) with OFF-DIAGONAL shifted means added:
D(a,b) = <sum_g Z(g+a)Z(g+b)>/<Z^2> = sinc(a-b) - sinc(a)sinc(b) (CUE limit, dab.py). Rows: m_0..m_J, d, D(al,al) (explore58's), D(a,b) pairs.
usage: cueAB.py N S delta J alphas pairs(a:b;a:b) [scales]
"""
import sys, math, json
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreLK')
from lk_cue import cue_angles, hardcore
PI = math.pi
def fvals(xs, P, pts):
    lf = np.sum(np.log(np.abs(np.sin(PI*(pts[:, None]-xs[None, :])/P))+1e-300), axis=1)
    sg = np.prod(np.sign(np.sin(PI*(pts[:, None]-xs[None, :])/P)), axis=1)
    return sg, lf
def column(xs, P, s, J, alphas, pairs):
    n2 = len(xs); n = n2//2; M = 4*n2+8
    t = np.arange(M)*(P/M)
    sg, lf = fvals(xs, P, t); sh = np.max(lf); f = sg*np.exp(lf-sh)
    c = np.fft.fft(f)/M; p = np.abs(c)**2
    k = np.rint(np.fft.fftfreq(M, d=1.0/M)).astype(int)
    w = np.zeros(n+1)
    for kk, pp in zip(k, p):
        if abs(kk) <= n: w[abs(kk)] += pp
    m0 = w.sum(); om = s*2*PI*np.arange(n+1)/P
    mj = [np.sum(w*om**(2*j))/m0 for j in range(J+1)]
    S = np.sin(PI*(xs[:, None]-xs[None, :])/P); np.fill_diagonal(S, 1.0)
    ls = np.sum(np.log(np.abs(S)), axis=1)-sh
    d = np.sum((PI/P)**2*np.exp(2*ls))/P/m0*s**3
    Ds = []
    for a in alphas:
        sgv, lv = fvals(xs, P, xs+s*a); Ds.append(s/P*np.sum(np.exp(2*(lv-sh)))/m0)
    Dp = []
    for (a, b) in pairs:
        sa, la = fvals(xs, P, xs+s*a); sb, lb = fvals(xs, P, xs+s*b)
        Dp.append(s/P*np.sum(sa*sb*np.exp(la+lb-2*sh))/m0)
    return np.array(mj+[d]+Ds+Dp)
def target(J, alphas, pairs):
    return np.array([PI**(2*j)/(2*j+1) for j in range(J+1)]+[PI**2/3]+[1-np.sinc(a)**2 for a in alphas]
                    +[np.sinc(a-b)-np.sinc(a)*np.sinc(b) for (a, b) in pairs])
def lp(V, t):
    A = V.T; R = 1.0/np.maximum(np.abs(t), 0.05); A = A*R[:, None]; b = t*R
    nr, nc = A.shape
    res = linprog(np.concatenate([np.zeros(nc), np.ones(2*nr)]), A_eq=np.hstack([A, np.eye(nr), -np.eye(nr)]), b_eq=b,
                  bounds=[(0, None)]*(nc+2*nr), method='highs')
    r = res.x[nc:nc+nr]-res.x[nc+nr:]
    return res.fun/nr, r, res.x[:nc], res
if __name__ == '__main__':
    N = int(sys.argv[1]); S = int(sys.argv[2]); delta = float(sys.argv[3]); J = int(sys.argv[4])
    alphas = [float(a) for a in sys.argv[5].split(',')] if sys.argv[5] != '-' else []
    pairs = [tuple(float(x) for x in pr.split(':')) for pr in sys.argv[6].split(';')] if sys.argv[6] != '-' else []
    scales = [float(x) for x in sys.argv[7].split(',')] if len(sys.argv) > 7 else [1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5]
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 97
    rng = np.random.default_rng(seed)
    t = target(J, alphas, pairs)
    cols = []; mg = 1e9
    for i in range(S):
        xs = cue_angles(N, rng)*N/(2*PI)
        ys, g = hardcore(xs, float(N), delta+1e-6) if delta > 0 else (xs, 0)
        mg = min(mg, g)
        for s in scales: cols.append(column(ys, float(N), s, J, alphas, pairs))
    cols.append(column(np.arange(N)+0.5, float(N), 1.0, J, alphas, pairs))
    V = np.array(cols)
    val, r, x, res = lp(V, t)
    print(f'N={N} S={S} delta={delta} J={J} alphas={alphas} pairs={pairs} scales={scales} min gap={mg:.4f}')
    print(f'  target={np.round(t,4).tolist()}')
    print(f'  mean|rel resid|={val:.3e} cols used={int(np.sum(x>1e-10))}; resid per row={np.round(r,5).tolist()}', flush=True)
    np.save(f'/Users/adriansager/playground/riemann/.scratch/ceilaudit7/V_N{N}_S{S}_d{delta}.npy', V)
