# Structure of the optimal fan for a saved set: cone labels of the seed (first k) points, angles between fan axes and
# seed directions, apex offset.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_survey import d2
from p6_cones import rot, SIMP
F = np.load(sys.argv[1]); k = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
n = len(F); D = np.empty((n, n), np.float32)
for i0 in range(0, n, 500): D[i0:i0 + 500] = np.sqrt(np.maximum(d2(F[i0:i0 + 500], F), 0))
def obj(z):
    Q = rot(z[:6]); c = z[6:]; lab = np.argmax((F - c) @ (SIMP @ Q.T).T, 1); m = 0.0
    for j in range(5):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = max(m, float(D[np.ix_(idx, idx)].max()))
    return m
best = (9, None)
for rs in range(5):
    z = np.r_[rng.uniform(-np.pi, np.pi, 6), np.zeros(4)]; f = obj(z); sig = np.r_[np.full(6, 0.3), np.full(4, 0.03)]
    for it in range(300):
        Z = z[None] + sig[None] * rng.normal(size=(10, 10)); fs = [obj(q) for q in Z]; o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.93
        if sig[0] < 1e-3: break
    if f < best[0]: best = (f, z)
f, z = best; Q = rot(z[:6]); c = z[6:]; Ax = SIMP @ Q.T
lab = np.argmax((F[:k] - c) @ Ax.T, 1)
ang = np.degrees(np.arccos(np.clip((F[:k] / np.linalg.norm(F[:k], axis=1, keepdims=True)) @ Ax.T, -1, 1)))
print(json.dumps({'best': round(f, 4), 'apex': np.round(c, 4).tolist(), 'seed_labels': lab.tolist(), 'seed_axis_angles_deg': np.round(ang, 1).tolist()}))
