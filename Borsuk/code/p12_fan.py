# Fan test: build a nearly complete diameter-1 set F (greedy extreme-point completion of a seed inside the candidate
# pool), then minimise over simplex fans (rotation Q, apex c) the max piece diameter of F.  Seeds: 'ty' = two orthogonal
# unit triangles (type-3 extreme, rho = 1/3), 'c3' / 'c4' / 'c5' = random admissible contact sets (all on the circumsphere).
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_survey import random_config, d2
from p6_cones import rot, SIMP
seed_kind, seed, nadd, iters = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
if seed_kind == 'ty':
    R = 1 / np.sqrt(3); th = 2 * np.pi * np.arange(3) / 3
    F = np.concatenate([R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1), R * np.stack([0 * th, 0 * th, np.cos(th + np.pi / 2), np.sin(th + np.pi / 2)], 1)])
    rho = 1 / 3
else:
    rho, P, cuts, lam = random_config(int(seed_kind[1]), rng); F = P.copy(); R = np.sqrt(rho)
# candidate pool inside B(0, R) within distance 1 of all seed points
pool = rng.normal(size=(400000, 4)); pool /= np.linalg.norm(pool, axis=1, keepdims=True); pool *= R * rng.uniform(0.3, 1, (len(pool), 1)) ** 0.25
ok = np.ones(len(pool), bool)
for f in F: ok &= ((pool - f) ** 2).sum(1) <= 1
pool = pool[ok]
t0 = time.time()
for k in range(nadd):
    if len(pool) == 0: break
    v = rng.normal(size=4); v /= np.linalg.norm(v); i = int(np.argmax(pool @ v)); z = pool[i]
    F = np.vstack([F, z]); pool = pool[((pool - z) ** 2).sum(1) <= 1]
n = len(F); D = np.empty((n, n), np.float32)
for i0 in range(0, n, 500): D[i0:i0 + 500] = np.sqrt(np.maximum(d2(F[i0:i0 + 500], F), 0))
print(json.dumps({'seed': seed_kind, 'rho': round(float(rho), 4), 'n': n, 'diam': round(float(D.max()), 5), 'pairs_ge_0.99': int(sum(int((D[i0:i0 + 500] >= 0.99).sum()) for i0 in range(0, n, 500)) // 2), 'pool_left': len(pool), 'sec': round(time.time() - t0)}), flush=True)
np.save('Borsuk/.scratch/fanset_%s_%d.npy' % (seed_kind, seed), F)
def obj(z):
    Q = rot(z[:6]); c = z[6:]; Dir = SIMP @ Q.T; lab = np.argmax((F - c) @ Dir.T, 1); m = 0.0
    for j in range(5):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = max(m, D[np.ix_(idx, idx)].max())
    return m
best = (9, None)
for rs in range(iters):
    z = np.r_[rng.uniform(-np.pi, np.pi, 6), 0.05 * rng.normal(size=4)]; f = obj(z); sig = np.r_[np.full(6, 0.3), np.full(4, 0.03)]
    for it in range(300):
        Z = z[None] + sig[None] * rng.normal(size=(12, 10)); fs = np.array([obj(q) for q in Z]); o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.93
        if sig[0] < 1e-3: break
    if f < best[0]: best = (f, z)
    print(json.dumps({'restart': rs, 'fan_max_piece': round(float(f), 5), 'best': round(float(best[0]), 5)}), flush=True)
