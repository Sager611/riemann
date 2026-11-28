# Adversarial fan game in R^4: X (diameter <= 1) starts from a seed; repeat: best simplex fan (rotation Q, apex c) for X;
# adversary adds a pool point z (within 1 of all of X) that has a partner x in X with |x - z| >= 1 - tol in the same cone
# of the current fan (choosing the z that keeps most pool points alive).  Track the best fan value per round.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_survey import random_config, d2
from p6_cones import rot, SIMP
seed_kind, seed, rounds, tol = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
rng = np.random.default_rng(seed)
if seed_kind == 'ty':
    R = 1 / np.sqrt(3); th = 2 * np.pi * np.arange(3) / 3
    X = np.concatenate([R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1), R * np.stack([0 * th, 0 * th, np.cos(th + np.pi / 2), np.sin(th + np.pi / 2)], 1)])
elif seed_kind == 'tri':
    R = 1 / np.sqrt(3); th = 2 * np.pi * np.arange(3) / 3; X = R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1)
else:
    rho, P, cuts, lam = random_config(int(seed_kind[1]), rng); X = P.copy(); R = np.sqrt(rho)
pool = rng.normal(size=(300000, 4)); pool /= np.linalg.norm(pool, axis=1, keepdims=True); pool *= R * rng.uniform(0.45, 1, (len(pool), 1)) ** 0.25
for x in X: pool = pool[((pool - x) ** 2).sum(1) <= 1]
for k in range(600):      # initial greedy completion (boundary points)
    v = rng.normal(size=4); z = pool[int(np.argmax(pool @ v))]; X = np.vstack([X, z]); pool = pool[((pool - z) ** 2).sum(1) <= 1]
def fanval(X, z):
    Q = rot(z[:6]); Dir = SIMP @ Q.T; lab = np.argmax((X - z[6:]) @ Dir.T, 1); m = 0.0; arg = None
    for j in range(5):
        Y = X[lab == j]
        if len(Y) > 1:
            M = d2(Y, Y); v = float(M.max())
            if v > m: m = v; arg = j
    return np.sqrt(m), lab, arg
def best_fan(X, z0s, n_iter=150):
    best = (9, None)
    for z in z0s:
        f = fanval(X, z)[0]; sig = np.r_[np.full(6, 0.2), np.full(4, 0.02)]
        for it in range(n_iter):
            Z = z[None] + sig[None] * rng.normal(size=(8, 10)); fs = [fanval(X, q)[0] for q in Z]; o = int(np.argmin(fs))
            if fs[o] < f: z, f = Z[o], fs[o]
            else: sig *= 0.9
            if sig[0] < 2e-3: break
        if f < best[0]: best = (f, z)
    return best
zbest = None; t0 = time.time()
for r in range(rounds):
    starts = [np.r_[rng.uniform(-np.pi, np.pi, 6), np.zeros(4)] for _ in range(3)] + ([zbest] if zbest is not None else [])
    f, zbest = best_fan(X, starts)
    # adversary: pool points z with a partner in X at distance >= 1 - tol in the same cone of the current fan
    Q = rot(zbest[:6]); Dir = SIMP @ Q.T; labX = np.argmax((X - zbest[6:]) @ Dir.T, 1); labP = np.argmax((pool - zbest[6:]) @ Dir.T, 1)
    cand = []
    for j in range(5):
        Pj = np.flatnonzero(labP == j); Xj = X[labX == j]
        if len(Pj) == 0 or len(Xj) == 0: continue
        for i0 in range(0, len(Pj), 2000):
            M = d2(pool[Pj[i0:i0 + 2000]], Xj).max(1); cand += Pj[i0:i0 + 2000][M >= (1 - tol) ** 2].tolist()
    print(json.dumps({'round': r, 'n': len(X), 'best_fan': round(float(f), 5), 'adversary_candidates': len(cand), 'pool': len(pool), 'sec': round(time.time() - t0)}), flush=True)
    if not cand: break
    cand = rng.choice(cand, min(60, len(cand)), replace=False)
    keep = [int((((pool - pool[c]) ** 2).sum(1) <= 1).sum()) for c in cand]; c = cand[int(np.argmax(keep))]
    z = pool[c]; X = np.vstack([X, z]); pool = pool[((pool - z) ** 2).sum(1) <= 1]
np.save('Borsuk/.scratch/adv_%s_%d.npy' % (seed_kind, seed), X)
