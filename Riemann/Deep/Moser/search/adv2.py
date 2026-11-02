import json, sys, time
import numpy as np
from scipy.optimize import minimize  # pyright: ignore[reportMissingImports]
from adversary import Shape, polyline
from shapeopt import to_params

cfg = json.load(open(sys.argv[1]))
S = Shape(cfg["t"], cfg["h"], ngrid=6000)
rng = np.random.default_rng(int(sys.argv[2]))
seeds = []
for f in sys.argv[3].split(','):
    if f:
        seeds.append(np.array(json.load(open(f))['Z']))
# U-shape family seeds with splayed legs
for a in np.linspace(0.18, 0.32, 8):
    for ang in np.radians(np.linspace(60, 120, 7)):
        s = 1 - 2 * a
        Z = np.array([[0, 0], [a, 0], [a + s * np.cos(ang), s * np.sin(ang)], [s * np.cos(ang), s * np.sin(ang)]])
        Z[3] = Z[2] + a * np.array([np.cos(np.pi - (np.pi - ang) * 0 - (np.pi - ang) - ang + np.pi), 0]) if False else Z[3]
        seeds.append(Z)
best = []
t0 = time.time()
def opt(x0, m, iters=3000):
    f = lambda x: -S.rho(polyline(x, m))
    res = minimize(f, x0, method='Nelder-Mead', options=dict(maxiter=iters, xatol=1e-9, fatol=1e-12, adaptive=True))
    return -res.fun, res.x
for Z in seeds:
    x0, m = to_params(Z)
    if m < 2:
        continue
    v, x = opt(x0, m)
    best.append((v, polyline(x, m)))
for m in (3, 4, 5, 6, 7, 8, 10):
    for _ in range(int(sys.argv[4])):
        x0 = np.concatenate([rng.uniform(-np.pi, np.pi, m - 1), rng.normal(0, 0.5, m)])
        v, x = opt(x0, m)
        best.append((v, polyline(x, m)))
best.sort(key=lambda t: -t[0])
S2 = Shape(cfg['t'], cfg['h'], ngrid=20000)
out = [{'rho': float(S2.rho(Z)), 'Z': Z.round(7).tolist()} for v, Z in best[:8]]
json.dump({'top': out, 'secs': round(time.time() - t0)}, open(sys.argv[5], 'w'), indent=1)
print(json.dumps([o['rho'] for o in out]))

