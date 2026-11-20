#!/usr/bin/env python3
"""explore58: DC-free real-rooted trig polynomials of type pi (2n zeros per period 2n).
Minimise e = d/m_1 (Gonek-type ratio; zeta: e = 1) subject to <f> = 0 (no atom at frequency 0) and all cyclic gaps >= delta.
"""
import sys, math, json
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, 'Riemann/Deep/ExploreLK')
from lk_cone import zeros_from_params
PI = math.pi

def spec(xs, P):
    n = len(xs) // 2; N = 8 * n + 16
    t = np.arange(N) * (P / N)
    f = np.prod(np.sin(PI * (t[:, None] - xs[None, :]) / P), axis=1)
    c = np.fft.fft(f) / N
    k = np.fft.fftfreq(N, d=1.0 / N)
    p = np.abs(c) ** 2
    w = 2 * PI * k / P
    m0 = p.sum(); m1 = np.sum(p * w ** 2); dc = p[0]
    S = np.sin(PI * (xs[:, None] - xs[None, :]) / P); np.fill_diagonal(S, 1.0)
    d = np.sum(((PI / P) * np.prod(S, axis=1)) ** 2) / P
    return m0, m1, dc, d, p, w

def best_e(n, delta, rng, restarts=30, pen=1e4):
    n2 = 2 * n; P = float(n2)
    if n2 * delta > P + 1e-12: return None
    def obj(z):
        xs = zeros_from_params(z, n2, P, delta)
        m0, m1, dc, d, p, w = spec(xs, P)
        return d / m1 + pen * dc / m0
    best = None
    for _ in range(restarts):
        z0 = rng.normal(size=n2) * rng.uniform(0.2, 4)
        r = minimize(obj, z0, method='Nelder-Mead', options={'maxiter': 3000 * n2, 'xatol': 1e-9, 'fatol': 1e-13})
        r = minimize(obj, r.x, method='Nelder-Mead', options={'maxiter': 3000 * n2, 'xatol': 1e-10, 'fatol': 1e-14})
        xs = zeros_from_params(r.x, n2, P, delta)
        m0, m1, dc, d, p, w = spec(xs, P)
        rec = dict(e=d / m1, dcfrac=dc / m0, m1=m1 / m0, d=d / m0, gaps=np.round(np.diff(np.concatenate([xs, [P]])), 5).tolist())
        if best is None or rec['e'] + pen * rec['dcfrac'] < best['e'] + pen * best['dcfrac']:
            best = rec
    return best

if __name__ == '__main__':
    rng = np.random.default_rng(7)
    out = {}
    for delta in [float(a) for a in sys.argv[1].split(',')]:
        for n in [int(a) for a in sys.argv[2].split(',')]:
            b = best_e(n, delta, rng, restarts=int(sys.argv[3]) if len(sys.argv) > 3 else 20)
            if b is None: continue
            print(f'delta={delta:.3f} n={n} min e={b["e"]:.6f} dc={b["dcfrac"]:.1e} m1/m0={b["m1"]:.4f} d/m0={b["d"]:.4f} gaps={b["gaps"]}', flush=True)
            out[f'{delta}_{n}'] = b
    json.dump(out, open('Riemann/Deep/ExploreLK/logs/dcfree_' + sys.argv[1].replace(',', '_') + '.json', 'w'), indent=1)
