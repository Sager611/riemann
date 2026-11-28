# 6-to-5 merging test.  A 6-piece partition {P_1..P_6} of an outer body Omega (pieces of diameter < 1) yields a 5-partition
# of every diameter-1 X ⊂ Omega unless X has unit pairs across all 15 piece pairs (conflict graph K_6).  Here: optimise a
# 6-cell power diagram on samples of Omega, then a greedy adversary tries to realise K_6 with unit pairs from a pool.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary, d2
specarg, seed, tol = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
spec = json.loads(open(specarg).read()) if specarg.endswith('.json') else json.loads(specarg)
rng = np.random.default_rng(seed); rho, P, cuts = body(spec); K = 6
S = np.concatenate([boundary(rho, P, cuts, 4000, rng), interior(rho, P, cuts, 1000, rng)])
DS = np.empty((len(S), len(S)), np.float32)
for i0 in range(0, len(S), 500): DS[i0:i0 + 500] = d2(S[i0:i0 + 500], S)
def lab_of(z, Y): G = z[:4 * K].reshape(K, 4); c = z[4 * K:]; return np.argmax(Y @ G.T + c[None], 1)
def phi(z):
    lab = lab_of(z, S); m = 0.0
    for j in range(K):
        idx = np.flatnonzero(lab == j)
        if len(idx) > 1: m = max(m, float(DS[np.ix_(idx, idx)].max()))
    return np.sqrt(m)
best = (9, None)
for rs in range(4):
    z = np.r_[rng.normal(size=4 * K), np.zeros(K)]; f = phi(z); sig = np.r_[np.full(4 * K, 0.3), np.full(K, 0.05)]
    for it in range(300):
        Z = z[None] + sig[None] * rng.normal(size=(10, 5 * K)); fs = [phi(q) for q in Z]; o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.95
        if sig[0] < 2e-3: break
    if f < best[0]: best = (f, z)
    print(json.dumps({'restart': rs, 'six_piece_max': round(f, 4), 'best': round(best[0], 4)}), flush=True)
z = best[1]
# adversary
del DS
pool = np.concatenate([boundary(rho, P, cuts, 40000, rng), interior(rho, P, cuts, 10000, rng)])
X = np.array(P, float); labX = list(lab_of(z, X)); conn = np.zeros((K, K), bool)
def update_conn(Xa, la):
    D = d2(Xa, Xa); a, b = np.nonzero(np.triu((D >= (1 - tol) ** 2) & (D <= 1 + 1e-12), 1))
    for i, j in zip(a, b):
        if la[i] != la[j]: conn[la[i], la[j]] = conn[la[j], la[i]] = True
for x in X: pool = pool[((pool - x) ** 2).sum(1) <= 1]
update_conn(X, labX); t0 = time.time()
for step in range(200):
    todo = [(i, j) for i in range(K) for j in range(i + 1, K) if not conn[i, j]]
    if not todo or len(pool) < 2: break
    lp = lab_of(z, pool); cands = []
    for i, j in todo:
        A = np.flatnonzero(lp == i); B = np.flatnonzero(lp == j)
        if len(A) == 0 or len(B) == 0: continue
        A = A if len(A) <= 1500 else rng.choice(A, 1500, replace=False); B = B if len(B) <= 1500 else rng.choice(B, 1500, replace=False)
        M = d2(pool[A], pool[B]); a, b = np.nonzero((M >= (1 - tol) ** 2) & (M <= 1.0))
        if len(a): s = rng.choice(len(a), min(20, len(a)), replace=False); cands += [(A[a[q]], B[b[q]]) for q in s]
    if not cands: break
    scores = [int(((((pool - pool[a]) ** 2).sum(1) <= 1) & (((pool - pool[b]) ** 2).sum(1) <= 1)).sum()) for a, b in cands]
    a, b = cands[int(np.argmax(scores))]; za, zb = pool[a].copy(), pool[b].copy()
    X = np.vstack([X, za, zb]); labX = list(lab_of(z, X)); pool = pool[(((pool - za) ** 2).sum(1) <= 1) & (((pool - zb) ** 2).sum(1) <= 1)]
    update_conn(X, labX)
    print(json.dumps({'step': step, 'connected_pairs': int(np.triu(conn, 1).sum()), 'of': 15, 'n': len(X), 'pool': len(pool), 'sec': round(time.time() - t0)}), flush=True)
print(json.dumps({'final_connected': int(np.triu(conn, 1).sum()), 'missing': [(i, j) for i in range(K) for j in range(i + 1, K) if not conn[i, j]], 'n': len(X)}), flush=True)
