# Set-cover adversary for the fan conjecture: grow X (diameter <= 1, inside B(0,R)) by unit pairs from a shell pool so that
# as many of K random fan rotations (apex 0) as possible contain a unit pair of X inside one cone.  Then optimise the fan
# (rotation + apex) for the final X.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_survey import random_config, d2
from p6_cones import rot, SIMP
kind, seed, K, tol = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
rng = np.random.default_rng(seed)
if kind == 'tri':
    R = 1 / np.sqrt(3); th = 2 * np.pi * np.arange(3) / 3; X = R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1)
else:
    rho, P, cuts, lam = random_config(int(kind[1]), rng); X = P.copy(); R = np.sqrt(rho)
rho = R * R; smin = (-(4 / 3) * R + np.sqrt((16 / 9) * rho - 4 * (rho - 1))) / 2
pool = rng.normal(size=(100000, 4)); pool /= np.linalg.norm(pool, axis=1, keepdims=True); pool *= rng.uniform(smin, R, (len(pool), 1))
for x in X: pool = pool[((pool - x) ** 2).sum(1) <= 1]
Qs = []
for k in range(K):
    A = rng.normal(size=(4, 4)); Q, Rr = np.linalg.qr(A); Qs.append(Q * np.sign(np.diag(Rr))[None, :])
Axes = np.stack([SIMP @ Q.T for Q in Qs])            # (K, 5, 4)
def labels(Y):
    out = np.empty((K, len(Y)), np.int8)
    for k in range(K): out[k] = np.argmax(Y @ Axes[k].T, 1)
    return out
def blocked(Y):
    L = labels(Y); D = d2(Y, Y); unit = (D >= (1 - tol) ** 2) & (D <= 1 + 1e-12); np.fill_diagonal(unit, False)
    b = np.zeros(K, bool)
    for k in range(K):
        same = L[k][:, None] == L[k][None, :]; b[k] = bool((unit & same).any())
    return b
t0 = time.time(); blk = blocked(X)
for step in range(400):
    if blk.all() or len(pool) < 2: break
    ks = rng.choice(np.flatnonzero(~blk), min(20, int((~blk).sum())), replace=False); cands = []
    Lp = labels(pool)
    for k in ks:
        for j in range(5):
            idx = np.flatnonzero(Lp[k] == j)
            if len(idx) < 2: continue
            sub = idx if len(idx) <= 1400 else rng.choice(idx, 1400, replace=False)
            M = d2(pool[sub], pool[sub]).astype(np.float32); a, b = np.nonzero(np.triu((M >= (1 - tol) ** 2) & (M <= 1.0), 1)); del M
            if len(a): s = rng.choice(len(a), min(15, len(a)), replace=False); cands += list(zip(sub[a[s]], sub[b[s]]))
    if not cands: break
    best = None
    for a, b in cands:
        la = Lp[:, a]; lb = Lp[:, b]; gain = int(((la == lb) & ~blk).sum())
        surv = int(((((pool - pool[a]) ** 2).sum(1) <= 1) & (((pool - pool[b]) ** 2).sum(1) <= 1)).sum())
        score = gain * np.sqrt(surv)
        if best is None or score > best[0]: best = (score, a, b, gain)
    _, a, b, gain = best; za, zb = pool[a].copy(), pool[b].copy()
    X = np.vstack([X, za, zb]); pool = pool[(((pool - za) ** 2).sum(1) <= 1) & (((pool - zb) ** 2).sum(1) <= 1)]
    blk |= (Lp[:, a] == Lp[:, b])
    if step % 10 == 0: print(json.dumps({'step': step, 'n': len(X), 'blocked': int(blk.sum()), 'of': K, 'pool': len(pool), 'sec': round(time.time() - t0)}), flush=True)
print(json.dumps({'final_n': len(X), 'blocked': int(blk.sum()), 'of': K, 'pool': len(pool), 'diam': round(float(np.sqrt(d2(X, X).max())), 5)}), flush=True)
np.save('Borsuk/.scratch/adv2_%s_%d.npy' % (kind, seed), X)
# free fan optimisation on the final X
D = np.sqrt(np.maximum(d2(X, X), 0))
def obj(z):
    Q = rot(z[:6]); c = z[6:]; lab = np.argmax((X - c) @ (SIMP @ Q.T).T, 1); m = 0.0
    for j in range(5):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = max(m, float(D[np.ix_(idx, idx)].max()))
    return m
bestf = 9
for rs in range(8):
    z = np.r_[rng.uniform(-np.pi, np.pi, 6), np.zeros(4)]; f = obj(z); sig = np.r_[np.full(6, 0.3), np.full(4, 0.03)]
    for it in range(300):
        Z = z[None] + sig[None] * rng.normal(size=(10, 10)); fs = [obj(q) for q in Z]; o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.93
        if sig[0] < 1e-3: break
    bestf = min(bestf, f); print(json.dumps({'restart': rs, 'fan': round(f, 5), 'best': round(bestf, 5)}), flush=True)
