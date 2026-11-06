#!/usr/bin/env python3
"""odlyzko31 optimiser (float): for each (n, r1) maximise [r1 logA + (n-r1) logB - E]/n over kernels
F = S(x/b) sech(wp x), S in {G,H}, b > 0, wp >= w. Writes a JSON of optimal kernels."""
import math, sys, json
import numpy as np
from kern import consts, EG, L8P

def grid(w, shapes=('G', 'H'), nb=260, nw=40, wmax_extra=1.6):
    bs = np.exp(np.linspace(math.log(0.05), math.log(400.0), nb))
    wps = w + np.concatenate([[0.0], np.geomspace(1e-3, wmax_extra, nw - 1)])
    rows = []
    for s in shapes:
        for b in bs:
            for wp in wps:
                lA, lB, E = consts(s, b, wp)
                if np.isfinite(E):
                    rows.append((s, b, wp, lA, lB, E))
    return rows

def best_on_grid(rows, n, r1):
    lA = np.array([r[3] for r in rows]); lB = np.array([r[4] for r in rows]); E = np.array([r[5] for r in rows])
    val = (r1 * lA + (n - r1) * lB - E) / n
    i = int(np.argmax(val))
    return rows[i], float(val[i])

def refine(n, r1, w, s, b0, wp0, iters=60):
    def f(lb, dw):
        wp = w + abs(dw)
        lA, lB, E = consts(s, math.exp(lb), wp)
        return (r1 * lA + (n - r1) * lB - E) / n
    # Nelder-Mead in (log b, dw), dw folded by abs
    pts = [np.array([math.log(b0), wp0 - w]), np.array([math.log(b0) + 0.05, wp0 - w]), np.array([math.log(b0), wp0 - w + 0.02])]
    vals = [f(*p) for p in pts]
    for _ in range(iters):
        order = np.argsort(vals)[::-1]; pts = [pts[i] for i in order]; vals = [vals[i] for i in order]
        c = (pts[0] + pts[1]) / 2
        xr = c + (c - pts[2]); fr = f(*xr)
        if fr > vals[0]:
            xe = c + 2 * (c - pts[2]); fe = f(*xe)
            pts[2], vals[2] = (xe, fe) if fe > fr else (xr, fr)
        elif fr > vals[1]:
            pts[2], vals[2] = xr, fr
        else:
            xc = c + 0.5 * (pts[2] - c); fc = f(*xc)
            if fc > vals[2]:
                pts[2], vals[2] = xc, fc
            else:
                pts[1] = pts[0] + 0.5 * (pts[1] - pts[0]); pts[2] = pts[0] + 0.5 * (pts[2] - pts[0])
                vals[1] = f(*pts[1]); vals[2] = f(*pts[2])
    i = int(np.argmax(vals))
    return math.exp(pts[i][0]), w + abs(pts[i][1]), vals[i]

if __name__ == '__main__':
    w = float(eval(sys.argv[1]))
    ns = [int(v) for v in sys.argv[2].split(',')]
    rows = grid(w)
    for n in ns:
        for r1 in ((0 if n % 2 == 0 else 1), n):
            (s, b, wp, *_), v = best_on_grid(rows, n, r1)
            b2, wp2, v2 = refine(n, r1, w, s, b, wp)
            print(f'n={n:6d} r1={r1:6d} shape={s} b={b2:9.4f} wp={wp2:.4f} (wp-w={wp2-w:.4f}) rd>{math.exp(v2):9.4f}  grid {math.exp(v):9.4f}')
