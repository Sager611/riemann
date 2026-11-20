#!/usr/bin/env python3
"""explore58: CUE characteristic polynomials as exact real-rooted type-pi trig polynomials (period N, density 1).
Target (discrete analogue of Z): flat spectrum over frequencies -N/2..N/2 and Gonek ratio d/m1 = 1 (zeta, RH).
Columns: CUE samples, optionally 'hard-core' modified to all gaps >= delta. LP feasibility (cone).
"""
import sys, math, json
import numpy as np
from scipy.optimize import linprog
PI = math.pi

def cue_angles(N, rng):
    Z = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R); Q = Q * (d / np.abs(d))
    ev = np.linalg.eigvals(Q)
    return np.sort(np.mod(np.angle(ev), 2 * PI))

def stats(xs, P):
    """f(x) = prod sin(pi (x - x_k)/P); returns folded spectrum weights w[0..n] (n = len/2) normalised to mass 1, d/m0, m1/m0."""
    n2 = len(xs); n = n2 // 2; N = 4 * n2 + 8
    t = np.arange(N) * (P / N)
    lf = np.sum(np.log(np.abs(np.sin(PI * (t[:, None] - xs[None, :]) / P))), axis=1)
    sg = np.prod(np.sign(np.sin(PI * (t[:, None] - xs[None, :]) / P)), axis=1)
    sh = np.max(lf); f = sg * np.exp(lf - sh)
    c = np.fft.fft(f) / N
    p = np.abs(c) ** 2
    k = np.rint(np.fft.fftfreq(N, d=1.0 / N)).astype(int)
    w = np.zeros(n + 1)
    for kk, pp in zip(k, p):
        if abs(kk) <= n: w[abs(kk)] += pp
    m0 = w.sum()
    S = np.sin(PI * (xs[:, None] - xs[None, :]) / P); np.fill_diagonal(S, 1.0)
    ls = np.sum(np.log(np.abs(S)), axis=1) - sh
    d = np.sum((PI / P) ** 2 * np.exp(2 * ls)) / P
    om = 2 * PI * np.arange(n + 1) / P
    return w / m0, d / m0, float(np.sum(w * om ** 2) / m0)

def hardcore(xs, P, delta, iters=2000):
    """push zeros apart until every cyclic gap >= delta (keeps order and period)."""
    x = xs.copy(); n = len(x)
    for _ in range(iters):
        g = np.diff(np.concatenate([x, [x[0] + P]]))
        bad = np.where(g < delta - 1e-12)[0]
        if len(bad) == 0: break
        for i in bad:
            j = (i + 1) % n
            gi = (x[j] - x[i]) % P
            sh = (delta - gi) / 2 + 1e-9
            x[i] -= sh; x[j] += sh
        x = np.sort(np.mod(x, P))
    g = np.diff(np.concatenate([x, [x[0] + P]]))
    return x, float(np.min(g))

def lp(cols, ds, tgt, dt):
    A = np.vstack([np.array(cols).T, np.array(ds)[None, :]])
    b = np.concatenate([tgt, [dt]])
    R = np.concatenate([1.0 / tgt, [1.0 / dt]]); A = A * R[:, None]; b = b * R
    nr, nc = A.shape
    c = np.concatenate([np.zeros(nc), np.ones(2 * nr)])
    res = linprog(c, A_eq=np.hstack([A, np.eye(nr), -np.eye(nr)]), b_eq=b, bounds=[(0, None)] * (nc + 2 * nr), method='highs')
    x = res.x[:nc]
    return res.fun / nr, int(np.sum(x > 1e-10)), res.x[nc:]

if __name__ == '__main__':
    N = int(sys.argv[1]); S = int(sys.argv[2]); deltas = [float(a) for a in sys.argv[3].split(',')]
    rng = np.random.default_rng(2026)
    n = N // 2; P = float(N)
    tgt = np.array([1.0] + [2.0] * n) / (N + 1)
    om = 2 * PI * np.arange(n + 1) / P
    m1t = float(np.sum(tgt * om ** 2)); dt = m1t  # Gonek ratio 1
    samples = [cue_angles(N, rng) * N / (2 * PI) for _ in range(S)]
    # natural CUE average (char-poly normalisation): flat spectrum check and d/m1
    W = []; D = []; G = []
    for xs in samples:
        w, d, m1 = stats(xs, P); W.append(w); D.append(d); G.append(np.min(np.diff(np.concatenate([xs, [xs[0] + P]]))))
    print(f'N={N} S={S} target m1={m1t:.4f} (pi^2/3={PI**2/3:.4f}); per-sample mean d/m1={np.mean(np.array(D)/np.array([np.sum(w*om**2) for w in W])):.4f}; median min gap={np.median(G):.4f}', flush=True)
    out = {}
    for delta in deltas:
        cols = []; ds = []; mg = []
        for xs in samples:
            ys, g = hardcore(xs, P, delta) if delta > 0 else (xs, None)
            w, d, m1 = stats(ys, P); cols.append(w); ds.append(d); mg.append(g if g is not None else 0)
        cols.append(np.eye(n + 1)[n]); ds.append(2 * PI ** 2)  # pure top sinusoid (gaps 1)
        val, nu, resid = lp(cols, ds, tgt, dt)
        print(f'delta={delta:.3f}: min gap after hard-core {np.min(mg):.4f}; LP mean|resid|={val:.3e} cols used={nu}; d-row resid={resid[-1]-resid[-1-1] if False else round(float(resid[n+1]-resid[2*(n+2)-1]),5)}; mean d/m1 of columns={np.mean(np.array(ds[:-1])/np.array([np.sum(c*om**2) for c in cols[:-1]])):.4f}', flush=True)
        out[str(delta)] = dict(resid=val, used=nu)
    json.dump(out, open(f'Riemann/Deep/ExploreLK/logs/cue_N{N}.json', 'w'), indent=1)
