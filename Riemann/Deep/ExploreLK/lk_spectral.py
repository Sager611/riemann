#!/usr/bin/env python3
"""explore58: all-moment version. Target: spectral measure of Z uniform on (0, pi] (no atom at 0), mass 1, and Gonek d = pi^2/3.
Columns: scaled copies (s in (0,1]) of DC-free real-rooted type-pi bases with all gaps >= delta (scaling only stretches gaps).
Binned spectral LP (K equal bins on (0, pi]) + d row; phase-1 infeasibility -> 0 means consistent.
"""
import sys, math, json
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, 'Riemann/Deep/ExploreLK')
from lk_dcfree import best_e, spec
PI = math.pi

def base_record(xs, P):
    m0, m1, dc, d, p, w = spec(xs, P)
    pos = w > 1e-12
    return dict(freqs=w[pos], wts=2 * p[pos] / m0, d=d / m0, dc=dc / m0)  # 2*: +-k pairs, normalised mass 1 on (0,pi]

def lp(bases, K, svals, dZ=PI ** 2 / 3):
    cols = []; meta = []
    for bi, b in enumerate(bases):
        for s in svals:
            col = np.zeros(K + 1)
            for f, wt in zip(b['freqs'], b['wts']):
                om = s * f
                k = min(int(om / PI * K), K - 1)
                col[k] += wt
            col[K] = b['d'] * s ** 3
            cols.append(col); meta.append((bi, float(s)))
    A = np.array(cols).T
    tgt = np.concatenate([np.full(K, 1.0 / K), [dZ]])
    R = np.concatenate([np.full(K, float(K)), [1.0 / dZ]])  # row scaling
    A = A * R[:, None]; tgt = tgt * R
    nr, nc = A.shape
    c = np.concatenate([np.zeros(nc), np.ones(2 * nr)])
    Aeq = np.hstack([A, np.eye(nr), -np.eye(nr)])
    res = linprog(c, A_eq=Aeq, b_eq=tgt, bounds=[(0, None)] * (nc + 2 * nr), method='highs')
    x = res.x[:nc]
    used = sorted([(float(x[i]), meta[i]) for i in range(nc) if x[i] > 1e-9], key=lambda u: -u[0])
    return res.fun / nr, used

if __name__ == '__main__':
    rng = np.random.default_rng(11)
    K = int(sys.argv[2]); out = {}
    for delta in [float(a) for a in sys.argv[1].split(',')]:
        bases = [dict(freqs=np.array([PI]), wts=np.array([1.0]), d=2 * PI ** 2, dc=0.0)]  # sine
        info = ['sine']
        for n in [int(a) for a in sys.argv[3].split(',')]:
            b = best_e(n, delta, rng, restarts=int(sys.argv[4]))
            if b is None: continue
            gaps = np.array(b['gaps']); xs = np.concatenate([[0.0], np.cumsum(gaps[:-1])])
            rec = base_record(xs, float(2 * n)); bases.append(rec); info.append((n, round(b['e'], 4), round(min(gaps), 4)))
        svals = np.linspace(0.002, 1.0, 500)
        val, used = lp(bases, K, svals)
        print(f'delta={delta:.3f} K={K} bases={info} mean|resid|={val:.3e} ncols_used={len(used)} top={[(round(w,4),m) for w,m in used[:5]]}', flush=True)
        out[str(delta)] = dict(resid=val, bases=info, used=used[:40])
    json.dump(out, open('Riemann/Deep/ExploreLK/logs/spectral_lp_' + sys.argv[1].replace(',', '_') + '.json', 'w'), indent=1, default=str)
