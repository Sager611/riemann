# P5: Tolmachev's 3-contact normal form H(u) (arXiv:2611.16775 Def. H, n = 4).  Lower bound for the
# 5-part Borsuk diameter of H(u): beta5 of a finite boundary sample (binary search, SAT).
import numpy as np, sys, time, gc, json
from pysat.solvers import Solver

def H_params(u):
    rho4, c4 = 0.4, 5 / 6
    t = u[0] * np.sqrt(rho4 - 1 / 3) - (1 - u[0]) * np.sqrt(1 / 12); rho = 1 / 3 + t * abs(t)
    h = (1 - u[1]) * np.sqrt((1 - c4) * rho) + u[1] * max(t, 0); tau = rho - h * h
    a0 = np.sqrt(max(0, 1 - 1 / (4 * tau))); a = (1 - u[2]) / 2 + u[2] * a0
    b = (1 - u[3]) * np.sqrt((1 - a) / 2) + u[3] * a
    st = np.sqrt(tau)
    P = np.array([[st, 0, h, 0], [st * (2 * a * a - 1), 2 * st * a * np.sqrt(1 - a * a), h, 0],
                  [st * (2 * b * b - 1), -2 * st * b * np.sqrt(1 - b * b), h, 0]])
    return rho, P

def boundary_sample(rho, P, N, rng):
    d = rng.normal(size=(N, 4)); d /= np.linalg.norm(d, axis=1, keepdims=True)
    s = np.full(N, np.sqrt(rho))
    for p in P:
        dp = d @ p; s = np.minimum(s, dp + np.sqrt(dp * dp + 1 - rho))
    pos = d[:, 3] > 0; s[pos] = np.minimum(s[pos], 0.5 / d[pos, 3])
    return d * s[:, None] * (1 - 1e-12)

def colourable(X, th, k=5):
    N = len(X); s = Solver(name='cd153')
    for v in range(N): s.add_clause([v * k + c + 1 for c in range(k)])
    first = None
    for i0 in range(0, N, 200):
        D2 = ((X[i0:i0 + 200, None, :] - X[None]) ** 2).sum(-1)
        ii, jj = np.nonzero(D2 >= th * th); m = jj > ii + i0; ii = ii[m] + i0; jj = jj[m]
        for a, b in zip(ii.tolist(), jj.tolist()):
            if first is None: first = a
            for c in range(k): s.add_clause([-(a * k + c + 1), -(b * k + c + 1)])
    if first is not None: s.add_clause([first * k + 1])
    res = s.solve(); s.delete(); return res

def beta5(X, lo=0.80, hi=1.10, steps=8):
    if not colourable(X, hi): return (hi, None)
    for _ in range(steps):
        mid = (lo + hi) / 2
        if colourable(X, mid): hi = mid
        else: lo = mid
    return (lo, hi)

if __name__ == '__main__':
    N = int(sys.argv[1]); grid = [float(x) for x in sys.argv[2].split(',')]; seed = int(sys.argv[3])
    rng = np.random.default_rng(seed)
    for u0 in grid:
        for u1 in grid:
            for u2 in grid:
                for u3 in grid:
                    u = (u0, u1, u2, u3); rho, P = H_params(u)
                    X = boundary_sample(rho, P, N, rng)
                    t0 = time.time(); lo, hi = beta5(X)
                    print(json.dumps({'u': u, 'rho': round(rho, 4), 'beta5': [round(lo, 4), hi and round(hi, 4)], 'sec': round(time.time() - t0, 1)}), flush=True)
