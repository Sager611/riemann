#!/usr/bin/env python3
"""explore58: closed-form Duffin–Schaeffer pair dictionary. g_{a,s}(x) = cos(s*pi*x) - a, 0<=a<1, 0<s<=1, plus constants.
Exact: m0 = a^2 + 1/2, m_j = (s pi)^{2j}/2 (j>=1), d = s^3 pi^2 (1 - a^2), min gap = 2 arccos(a)/(s pi).
LP: is v_Z(J) = (pi^{2j}/(2j+1))_{j<=J}, d = pi^2/3 in the cone of columns with min gap >= delta?  -> delta_0(J) for this dictionary.
"""
import sys, math, json
import numpy as np
from scipy.optimize import linprog
PI = math.pi

def cols(J, delta, na=60, ns=300):
    C = [np.concatenate([[1.0], np.zeros(J), [0.0]])]; M = [('const',)]
    for s in np.linspace(1e-3, 1.0, ns):
        amax = math.cos(min(PI, delta * s * PI / 2))  # need 2 arccos(a)/(s pi) >= delta  <=> a <= cos(delta s pi/2)
        for a in np.linspace(0.0, amax, na):
            v = [a * a + 0.5] + [(s * PI) ** (2 * j) / 2 for j in range(1, J + 1)] + [s ** 3 * PI ** 2 * (1 - a * a)]
            v = np.array(v); C.append(v / v[0]); M.append(('pair', float(a), float(s)))
    return np.array(C), M

def feas(J, delta):
    V, M = cols(J, delta)
    t = np.array([PI ** (2 * j) / (2 * j + 1) for j in range(J + 1)] + [PI ** 2 / 3])
    R = np.array([PI ** (-2 * j) * (2 * j + 1) for j in range(J + 1)] + [3 / PI ** 2])
    A = (V * R[None, :]).T; b = t * R
    nr, nc = A.shape
    c = np.concatenate([np.zeros(nc), np.ones(2 * nr)])
    res = linprog(c, A_eq=np.hstack([A, np.eye(nr), -np.eye(nr)]), b_eq=b, bounds=[(0, None)] * (nc + 2 * nr), method='highs')
    x = res.x[:nc]
    used = sorted([(float(x[i]), M[i]) for i in range(nc) if x[i] > 1e-10], key=lambda u: -u[0])
    return res.fun, used

if __name__ == '__main__':
    out = {}
    for J in [int(a) for a in sys.argv[1].split(',')]:
        lo, hi = 0.0, 1.0
        r0, _ = feas(J, 0.0)
        if r0 > 1e-9:
            print(f'J={J}: infeasible already at delta=0 (resid {r0:.2e})', flush=True); out[J] = None; continue
        for _ in range(14):
            mid = (lo + hi) / 2
            r, _ = feas(J, mid)
            if r < 1e-9: lo = mid
            else: hi = mid
        r, used = feas(J, lo)
        print(f'J={J}: delta_0(pair dictionary) in [{lo:.4f}, {hi:.4f}]; mixture at {lo:.4f}: {[(round(w,5), m) for w, m in used[:6]]}', flush=True)
        out[J] = dict(lo=lo, hi=hi, used=used[:10])
    json.dump(out, open('Riemann/Deep/ExploreLK/logs/pairs_lp.json', 'w'), indent=1, default=str)
