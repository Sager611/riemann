# Dimension-3 control for the fan conjecture: nearly complete diameter-1 sets in R^3 (greedy completion of a seed)
# and tetrahedral fans (rotation Q in SO(3), apex c) minimising the max piece diameter.
import numpy as np, sys, json, time
from scipy.linalg import expm
seed_kind, seed, nadd, iters = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * A @ B.T
TET = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
def rot(w):
    A = np.array([[0, w[0], w[1]], [-w[0], 0, w[2]], [-w[1], -w[2], 0]]); return expm(A)
if seed_kind == 'tri':
    R = 1 / np.sqrt(3); th = 2 * np.pi * np.arange(3) / 3; F = R * np.stack([np.cos(th), np.sin(th), 0 * th], 1)
elif seed_kind == 'tet':
    R = np.sqrt(3 / 8); F = R * TET
elif seed_kind == 'tri2':        # triangle + pole: x = (0, 0, R) is within 1 of the triangle (|x - p|^2 = 2/3)
    R = 1 / np.sqrt(3); th = 2 * np.pi * np.arange(3) / 3; F = np.vstack([R * np.stack([np.cos(th), np.sin(th), 0 * th], 1), [0, 0, R]])
else:   # random 4-point contact set on a sphere with 0 inside, pairwise <= 1
    while True:
        rho = rng.uniform(0.32, 0.375); R = np.sqrt(rho); v = rng.normal(size=(4, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True); F = R * v
        if d2(F, F).max() <= 1 + 1e-12:
            A = np.vstack([F.T, np.ones(4)]); lam = np.linalg.lstsq(A, np.r_[0, 0, 0, 1.0], rcond=None)[0]
            if np.all(lam > 0.03): break
pool = rng.normal(size=(300000, 3)); pool /= np.linalg.norm(pool, axis=1, keepdims=True); pool *= R * rng.uniform(0.2, 1, (len(pool), 1)) ** (1 / 3)
ok = np.ones(len(pool), bool)
for f in F: ok &= ((pool - f) ** 2).sum(1) <= 1
pool = pool[ok]
for k in range(nadd):
    if len(pool) == 0: break
    v = rng.normal(size=3); v /= np.linalg.norm(v); z = pool[int(np.argmax(pool @ v))]
    F = np.vstack([F, z]); pool = pool[((pool - z) ** 2).sum(1) <= 1]
n = len(F); D = np.empty((n, n), np.float32)
for i0 in range(0, n, 500): D[i0:i0 + 500] = np.sqrt(np.maximum(d2(F[i0:i0 + 500], F), 0))
print(json.dumps({'dim': 3, 'seed': seed_kind, 'rho': round(float(R * R), 4), 'n': n, 'diam': round(float(D.max()), 5), 'pool_left': len(pool)}), flush=True)
def obj(z):
    Q = rot(z[:3]); c = z[3:]; Dir = TET @ Q.T; lab = np.argmax((F - c) @ Dir.T, 1); m = 0.0
    for j in range(4):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = max(m, float(D[np.ix_(idx, idx)].max()))
    return m
best = 9
for rs in range(iters):
    z = np.r_[rng.uniform(-np.pi, np.pi, 3), 0.05 * rng.normal(size=3)]; f = obj(z); sig = np.r_[np.full(3, 0.3), np.full(3, 0.03)]
    for it in range(300):
        Z = z[None] + sig[None] * rng.normal(size=(10, 6)); fs = np.array([obj(q) for q in Z]); o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.93
        if sig[0] < 1e-3: break
    best = min(best, f)
    print(json.dumps({'restart': rs, 'fan': round(float(f), 5), 'best': round(float(best), 5)}), flush=True)
