#!/usr/bin/env python3
"""explore58: cone test E(delta; I) for real-rooted trig polynomials of type <= pi (normalised units).
Column = base f(x) = prod_{k=1}^{2n} sin(pi (x - x_k)/P), P = 2n (type pi, density 1), scaled x -> s x (type s*pi).
v(f) = (m_0..m_J, d), m_j = <f^(j)^2>, d = zero-sum of f'^2 per unit length. Target v_Z = (pi^{2j}/(2j+1), pi^2/3).
Phase-1 LP + column generation (pricing over zero positions with cyclic gaps >= delta*s).
"""
import sys, json, math
import numpy as np
from scipy.optimize import linprog, minimize

PI = math.pi

def base_moments(xs, P, J):
    """xs: zeros in [0,P) sorted, 2n of them. returns m (J+1,), d (per unit length), normalised later."""
    n2 = len(xs); n = n2 // 2
    N = 8 * n + 16
    t = np.arange(N) * (P / N)
    f = np.prod(np.sin(PI * (t[:, None] - xs[None, :]) / P), axis=1)
    c = np.fft.fft(f) / N
    k = np.fft.fftfreq(N, d=1.0 / N)  # integer indices
    w = 2 * PI * k / P
    p = np.abs(c) ** 2
    m = np.array([np.sum(p * w ** (2 * j)) for j in range(J + 1)])
    # slopes at zeros
    diff = PI * (xs[:, None] - xs[None, :]) / P
    S = np.sin(diff)
    np.fill_diagonal(S, 1.0)
    slopes = (PI / P) * np.prod(S, axis=1)
    d = np.sum(slopes ** 2) / P
    return m, d

def column(xs, P, s, J):
    m, d = base_moments(xs, P, J)
    mj = m * s ** (2 * np.arange(J + 1))
    dd = d * s ** 3
    v = np.concatenate([mj, [dd]])
    return v / v[0]

def target(J):
    return np.array([PI ** (2 * j) / (2 * j + 1) for j in range(J + 1)] + [PI ** 2 / 3])

def rowscale(J):
    return np.array([PI ** (-2 * j) for j in range(J + 1)] + [PI ** -2])

def gaps_ok(xs, P, minbase):
    g = np.diff(np.concatenate([xs, [xs[0] + P]]))
    return np.min(g) >= minbase - 1e-12

def zeros_from_params(z, n2, P, gmin):
    # cyclic gaps >= gmin, sum = P; first zero at 0
    free = P - n2 * gmin
    e = np.exp(z - np.max(z)); e = e / e.sum()
    g = gmin + free * e
    xs = np.concatenate([[0.0], np.cumsum(g[:-1])])
    return xs

def phase1(V, t, R):
    A = (V * R[None, :]).T  # rows = moments
    b = t * R
    nc = A.shape[1]; nr = A.shape[0]
    c = np.concatenate([np.zeros(nc), np.ones(2 * nr)])
    Aeq = np.hstack([A, np.eye(nr), -np.eye(nr)])
    res = linprog(c, A_eq=Aeq, b_eq=b, bounds=[(0, None)] * (nc + 2 * nr), method='highs')
    y = res.eqlin.marginals if res.status == 0 else None
    return res.fun, res.x[:nc], y

def price(y, R, delta, J, rng, nmax=6, tries=40):
    """maximise y.(R*v(f)) over admissible columns; returns best list of (score, xs, P, s)."""
    best = []
    for n in range(1, nmax + 1):
        n2 = 2 * n; P = float(n2)
        for s in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]:
            gmin = delta * s
            if n2 * gmin > P: continue
            def obj(z):
                xs = zeros_from_params(z, n2, P, gmin)
                v = column(xs, P, s, J)
                return -float(np.dot(y, v * R))
            for _ in range(max(2, tries // nmax)):
                z0 = rng.normal(size=n2) * rng.uniform(0.1, 3)
                r = minimize(obj, z0, method='Nelder-Mead', options={'maxiter': 400 * n2, 'xatol': 1e-6, 'fatol': 1e-12})
                xs = zeros_from_params(r.x, n2, P, gmin)
                best.append((-r.fun, xs, P, s))
    best.sort(key=lambda b: -b[0])
    return best[:12]

def run(delta, J, seed=1, iters=25, verbose=False):
    rng = np.random.default_rng(seed)
    t = target(J); R = rowscale(J)
    cols = [np.concatenate([[1.0], np.zeros(J), [0.0]])]  # constant
    meta = [('const',)]
    # seed: pairs and sinusoids at several s
    for s in np.linspace(0.05, 1.0, 40):
        for a_g in np.linspace(delta * s, 1.0, 12):
            if a_g < delta * s - 1e-12: continue
            xs = np.array([0.0, a_g]); P = 2.0  # pair with gaps a_g, 2-a_g
            if 2 - a_g < delta * s: continue
            cols.append(column(xs, P, s, J)); meta.append(('pair', float(a_g), float(s)))
    hist = []
    for it in range(iters):
        V = np.array(cols)
        val, w, y = phase1(V, t, R)
        hist.append(val)
        if verbose: print(f'  it {it} infeas {val:.3e} ncols {len(cols)}', flush=True)
        if val < 1e-10: break
        cand = price(y, R, delta, J, rng)
        added = 0
        for sc, xs, P, s in cand:
            if sc > 1e-9:
                cols.append(column(xs, P, s, J)); meta.append(('gen', xs.tolist(), P, s)); added += 1
        if added == 0: break
    V = np.array(cols)
    val, w, y = phase1(V, t, R)
    used = [(float(w[i]), meta[i]) for i in range(len(w)) if w[i] > 1e-12]
    return val, used, y

if __name__ == '__main__':
    J = int(sys.argv[1]); deltas = [float(a) for a in sys.argv[2].split(',')]
    out = {}
    for dl in deltas:
        val, used, y = run(dl, J, verbose=False)
        print(f'J={J} delta={dl:.4f} infeas={val:.3e} used={len(used)} dual={np.round(y,4).tolist() if y is not None else None}', flush=True)
        out[str(dl)] = {'infeas': val, 'used': used[:10]}
    with open(f'Riemann/Deep/ExploreLK/logs/cone_J{J}.json', 'w') as fh:
        json.dump(out, fh, indent=1, default=str)
