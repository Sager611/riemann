import json, sys, time
import numpy as np
from scipy.optimize import minimize  # pyright: ignore[reportMissingImports]
from adversary import Shape, polyline
from shapeopt import to_params

cfg = json.load(open(sys.argv[1]))
S = Shape(cfg['t'], cfg['h'], ngrid=6000)
S2 = Shape(cfg['t'], cfg['h'], ngrid=40000)
top = json.load(open(sys.argv[2]))['top']
rng = np.random.default_rng(5)
res_all = []
t0 = time.time()
for o in top[:4]:
    Z = np.array(o['Z'])
    for sub in (2, 3):
        # subdivide each segment into 'sub' pieces and perturb, then re-optimize (lets legs bend)
        pts = [Z[0]]
        for a, b in zip(Z[:-1], Z[1:]):
            for s in range(1, sub + 1):
                pts.append(a + (b - a) * s / sub)
        W = np.array(pts)
        for trial in range(3):
            x0, m = to_params(W)
            x0 = x0 + rng.normal(0, 0.02, len(x0))
            f = lambda x: -S.rho(polyline(x, m))
            r = minimize(f, x0, method='Nelder-Mead', options=dict(maxiter=6000, xatol=1e-10, fatol=1e-13, adaptive=True))
            Zn = polyline(r.x, m)
            res_all.append((float(S2.rho(Zn)), Zn))
res_all.sort(key=lambda t: -t[0])
json.dump({'top': [{'rho': v, 'Z': Z.round(7).tolist()} for v, Z in res_all[:6]], 'secs': round(time.time() - t0)}, open(sys.argv[3], 'w'), indent=1)
print([round(v, 6) for v, _ in res_all[:6]])

