# Spherical set-cover adversary: grow D ⊂ S^3 (angular diameter <= theta, seeded so that 0 ∈ conv D) by pairs at angle
# in [theta - tol, theta] lying in one cone, to block as many of K random fan rotations as possible.  Then optimise the fan.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_cones import rot, SIMP
theta, seed, K, tol, kind = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
rng = np.random.default_rng(seed); c = np.cos(np.radians(theta)); c2 = np.cos(np.radians(theta - tol))
if kind == 'tri':
    t = 2 * np.pi * np.arange(3) / 3; D = np.stack([np.cos(t), np.sin(t), 0 * t, 0 * t], 1)
else:
    while True:
        v = rng.normal(size=(5, 4)); v /= np.linalg.norm(v, axis=1, keepdims=True)
        if (v @ v.T)[np.triu_indices(5, 1)].min() < c: continue
        A = np.vstack([v.T, np.ones(5)]); lam = np.linalg.lstsq(A, np.r_[0, 0, 0, 0, 1.0], rcond=None)[0]
        if np.allclose(A @ lam, np.r_[0, 0, 0, 0, 1.0]) and lam.min() > 0.02: D = v; break
pool = rng.normal(size=(100000, 4)); pool /= np.linalg.norm(pool, axis=1, keepdims=True)
for d in D: pool = pool[pool @ d >= c]
Ax = []
for k in range(K):
    A = rng.normal(size=(4, 4)); Q, Rr = np.linalg.qr(A); Ax.append(SIMP @ (Q * np.sign(np.diag(Rr))[None, :]).T)
def labels(Y):
    out = np.empty((K, len(Y)), np.int8)
    for k in range(K): out[k] = np.argmax(Y @ Ax[k].T, 1)
    return out
def blocked(Y):
    L = labels(Y); G = Y @ Y.T; pr = (G <= c2) & (G >= c - 1e-12); np.fill_diagonal(pr, False); b = np.zeros(K, bool)
    for k in range(K): b[k] = bool((pr & (L[k][:, None] == L[k][None, :])).any())
    return b
blk = blocked(D); t0 = time.time()
for step in range(500):
    if blk.all() or len(pool) < 2: break
    Lp = labels(pool); ks = rng.choice(np.flatnonzero(~blk), min(20, int((~blk).sum())), replace=False); cands = []
    for k in ks:
        for j in range(5):
            idx = np.flatnonzero(Lp[k] == j)
            if len(idx) < 2: continue
            sub = idx if len(idx) <= 1500 else rng.choice(idx, 1500, replace=False)
            G = (pool[sub] @ pool[sub].T).astype(np.float32); a, b = np.nonzero(np.triu((G <= c2) & (G >= c), 1)); del G
            if len(a): s = rng.choice(len(a), min(15, len(a)), replace=False); cands += list(zip(sub[a[s]], sub[b[s]]))
    if not cands: break
    best = None
    for a, b in cands:
        gain = int(((Lp[:, a] == Lp[:, b]) & ~blk).sum()); surv = int(((pool @ pool[a] >= c) & (pool @ pool[b] >= c)).sum())
        sc = gain * np.sqrt(surv)
        if best is None or sc > best[0]: best = (sc, a, b)
    _, a, b = best; za, zb = pool[a].copy(), pool[b].copy(); blk |= (Lp[:, a] == Lp[:, b])
    D = np.vstack([D, za, zb]); pool = pool[(pool @ za >= c) & (pool @ zb >= c)]
    if step % 10 == 0: print(json.dumps({'step': step, 'n': len(D), 'blocked': int(blk.sum()), 'of': K, 'pool': len(pool), 'sec': round(time.time() - t0)}), flush=True)
print(json.dumps({'final_n': len(D), 'blocked': int(blk.sum()), 'of': K, 'pool': len(pool)}), flush=True)
G = D @ D.T
def maxang(Q):
    lab = np.argmax(D @ (SIMP @ Q.T).T, 1); m = 1.0
    for j in range(5):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = min(m, float(G[np.ix_(idx, idx)].min()))
    return np.degrees(np.arccos(np.clip(m, -1, 1)))
bestf = 999
for rs in range(6):
    z = rng.uniform(-np.pi, np.pi, 6); f = maxang(rot(z)); sig = 0.3
    for it in range(250):
        Z = z[None] + sig * rng.normal(size=(8, 6)); fs = [maxang(rot(q)) for q in Z]; o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.92
        if sig < 1e-3: break
    bestf = min(bestf, f); print(json.dumps({'restart': rs, 'fan_maxang': round(f, 3), 'best': round(bestf, 3), 'theta': theta}), flush=True)
np.save('Borsuk/.scratch/sphadv_%s_%d.npy' % (kind, seed), D)
