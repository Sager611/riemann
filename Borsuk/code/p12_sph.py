# Spherical fan problem: D ⊂ S^3 with angular diameter <= theta and 0 in conv(D) (seeded with a contact configuration),
# greedily completed on the sphere.  Fraction of random fan rotations with no pair of D at angle >= theta - tol inside a cone;
# best fan by optimisation (max angle inside cones).
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_cones import rot, SIMP
theta = np.radians(float(sys.argv[1])); seed = int(sys.argv[2]); nadd = int(sys.argv[3]); kind = sys.argv[4]
rng = np.random.default_rng(seed); c = np.cos(theta)
def seedset():
    if kind == 'tri':   # great-circle triangle at 120 deg (only admissible if theta >= 120)
        t = 2 * np.pi * np.arange(3) / 3; return np.stack([np.cos(t), np.sin(t), 0 * t, 0 * t], 1)
    for _ in range(100000):   # random 5-point set with pairwise angle <= theta, 0 inside conv
        v = rng.normal(size=(5, 4)); v /= np.linalg.norm(v, axis=1, keepdims=True)
        if (v @ v.T)[np.triu_indices(5, 1)].min() < c: continue
        A = np.vstack([v.T, np.ones(5)]); lam = np.linalg.lstsq(A, np.r_[0, 0, 0, 0, 1.0], rcond=None)[0]
        if np.allclose(A @ lam, np.r_[0, 0, 0, 0, 1.0]) and lam.min() > 0.02: return v
D = seedset(); pool = rng.normal(size=(200000, 4)); pool /= np.linalg.norm(pool, axis=1, keepdims=True)
for d in D: pool = pool[pool @ d >= c]
for k in range(nadd):
    if len(pool) == 0: break
    v = rng.normal(size=4); z = pool[int(np.argmax(pool @ v))]; D = np.vstack([D, z]); pool = pool[pool @ z >= c]
G = D @ D.T; n = len(D)
def maxang(Q):
    lab = np.argmax(D @ (SIMP @ Q.T).T, 1); m = 1.0
    for j in range(5):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = min(m, float(G[np.ix_(idx, idx)].min()))
    return np.degrees(np.arccos(np.clip(m, -1, 1)))
vals = []
for k in range(400):
    A = rng.normal(size=(4, 4)); Q, Rr = np.linalg.qr(A); Q = Q * np.sign(np.diag(Rr))[None, :]; vals.append(maxang(Q))
vals = np.array(vals); th = np.degrees(theta)
best = 999
for rs in range(4):
    z = rng.uniform(-np.pi, np.pi, 6); f = maxang(rot(z)); sig = 0.3
    for it in range(250):
        Z = z[None] + sig * rng.normal(size=(8, 6)); fs = [maxang(rot(q)) for q in Z]; o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.92
        if sig < 1e-3: break
    best = min(best, f)
print(json.dumps({'theta': th, 'kind': kind, 'n': n, 'diam_deg': round(float(np.degrees(np.arccos(G.min()))), 3), 'good_frac(<theta-0.05)': round(float((vals < th - 0.05).mean()), 4),
                  'good_frac(<theta-1)': round(float((vals < th - 1).mean()), 4), 'median_maxang': round(float(np.median(vals)), 2), 'best_fan_maxang': round(float(best), 3)}))
