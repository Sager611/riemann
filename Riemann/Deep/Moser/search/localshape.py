import json, sys, time
import numpy as np
from scipy.optimize import minimize  # pyright: ignore[reportMissingImports]
from shapeopt import make_shape, poly_area, to_params, canonical
from adversary import polyline

# family: top arcs found for K0 (hard U's, Z's, hooks) + canonical
fam = []
for f, key in (('adv3_K0.json', 'top'), ('adv2b_K0.json', 'top'), ('adv2_K0.json', 'top')):
    for o in json.load(open(f))[key][:6]:
        fam.append(np.array(o['Z']))
for f in ['adv_K0_m3.json', 'adv_K0_m4.json', 'adv_K0_m5.json', 'adv_K0_m6.json']:
    fam.append(np.array(json.load(open(f))['Z']))
fam += canonical()
p0 = np.array([190.5034657088439, 47.25875546131363, 66.90661890814535, 0.852937690678452])


def F(p, fam):
    S = make_shape(p, 4, ngrid=6000)
    A = poly_area(S)
    r = np.array([S.rho(Z) for Z in fam])
    return A * r.max() ** 2, r


f0, r0 = F(p0, fam)
order = np.argsort(-r0)
print('F(K0)', f0, 'top rhos (normalized by h2=1)', (r0[order[:6]]).round(6).tolist())
# gradient by central differences
steps = np.array([0.05, 0.05, 0.05, 0.001])
g = []
for i in range(4):
    e = np.zeros(4); e[i] = steps[i]
    fp, _ = F(p0 + e, fam); fm, _ = F(p0 - e, fam)
    g.append((fp - fm) / (2 * steps[i]))
print('grad', np.array(g).tolist())
# NM on fixed family, then re-optimize the top arcs at the new shape
res = minimize(lambda p: F(p, fam)[0], p0, method='Nelder-Mead', options=dict(maxiter=300, initial_simplex=[p0] + [p0 + np.eye(4)[i] * steps[i] * 10 for i in range(4)], xatol=1e-6, fatol=1e-12))
p1 = res.x
f1, r1 = F(p1, fam)
print('NM fixed-family F', f1, 'p1', p1.tolist(), 'rel gain', 1 - f1 / f0)
S1 = make_shape(p1, 4, ngrid=6000)
fam2 = []
for Z in fam:
    x0, m = to_params(Z)
    if m < 2:
        fam2.append(Z); continue
    r = minimize(lambda x: -S1.rho(polyline(x, m)), x0, method='Nelder-Mead', options=dict(maxiter=1500, xatol=1e-9, fatol=1e-12, adaptive=True))
    Zn = polyline(r.x, m)
    fam2.append(Zn if S1.rho(Zn) > S1.rho(Z) else Z)
f1b, r1b = F(p1, fam2)
print('after re-opt of arcs at p1: F', f1b, 'rel gain vs K0', 1 - f1b / f0)
json.dump({'p1': p1.tolist(), 'F1': f1b, 'F0': f0, 'family': [Z.tolist() for Z in fam2]}, open('localshape.json', 'w'))

