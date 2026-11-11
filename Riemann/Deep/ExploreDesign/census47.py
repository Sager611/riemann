"""Census of the c = 47 forced lattice nodes beyond K2 (WeilZeta full rigorous scan, logs/scan_zeta_47_3683_rig.json)."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from wsym import mu, sigma
c = 47; L = math.log(c)
d = json.load(open('../WeilZeta/logs/scan_zeta_47_3683_rig.json'))
ko = np.array(d['selo'], dtype=float); je = np.array(d['sele'], dtype=float)
to = math.pi * ko / L; te = 2 * math.pi * je / L
mo = mu(to, c); me = mu(te, c)
print('sigma', sigma(c), 'odd', len(ko), 'even', len(je))
print('odd mu range', mo.min(), mo.max(), ' even mu range', me.min(), me.max())
edges = [3e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10]
for a, b in zip(edges[:-1], edges[1:]):
    so = (to >= a) & (to < b); se = (te >= a) & (te < b)
    print(f't in [{a:.0e},{b:.0e}): odd {so.sum():5d} (min mu {mo[so].min() if so.any() else float("nan"):6.2f})  even {se.sum():5d} (min mu {me[se].min() if se.any() else float("nan"):6.2f})  log(t/2pi) {math.log(a/2/math.pi):5.2f}')
# clusters: merge all forced t (odd and even) within gap 3
allt = np.sort(np.concatenate([to, te]))
gaps = np.diff(allt)
nclust = 1 + int((gaps > 3.0).sum())
print('clusters (gap>3):', nclust, ' median gap', float(np.median(gaps)))
np.save('/tmp/e27_forced47.npy', np.stack([np.concatenate([to, te]), np.concatenate([mo, me]), np.concatenate([np.ones_like(to), np.zeros_like(te)])]))
