# FC-k probe: k random fan rotations (apex 0) at rho; for each (rotation i, cone e) sample "zone pairs" (x, y in the cone,
# inside B(0,R), |x - y| in [1, 1+eps]); types (i,e), (j,f) are compatible if some zone pair of each has all 4 cross
# distances <= 1 + eps.  CSP: pick one cone per rotation, pairwise compatible (pysat).  UNSAT = no X blocks all k (sampled).
import numpy as np, sys, json, time
from pysat.solvers import Solver
sys.path.insert(0, 'Borsuk/code')
from p6_cones import SIMP
rho, k, npair, eps, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
R = np.sqrt(rho); rng = np.random.default_rng(seed)
smin = (-(4 / 3) * R + np.sqrt((16 / 9) * rho - 4 * (rho - 1))) / 2
def shell(n):
    v = rng.normal(size=(n, 4)); v /= np.linalg.norm(v, axis=1, keepdims=True); return v * ((rng.uniform(smin ** 4, R ** 4, n)) ** 0.25)[:, None]
Qs = []
for i in range(k):
    A = rng.normal(size=(4, 4)); Q, Rr = np.linalg.qr(A); Qs.append(Q * np.sign(np.diag(Rr))[None, :])
S = shell(60000); types = []; pairs = []
t0 = time.time()
for i, Q in enumerate(Qs):
    lab = np.argmax(S @ (SIMP @ Q.T).T, 1)
    for e in range(5):
        idx = np.flatnonzero(lab == e); P = S[idx]; got = []
        for i0 in range(0, len(P), 500):
            A = P[i0:i0 + 500]; D = (A * A).sum(1)[:, None] + (P * P).sum(1)[None] - 2 * A @ P.T
            a, b = np.nonzero((D >= 1) & (D <= (1 + eps) ** 2)); got += list(zip((i0 + a).tolist(), b.tolist()))
            if len(got) > 20 * npair: break
        if got:
            sel = rng.choice(len(got), min(npair, len(got)), replace=False); pr = np.array([[P[got[s][0]], P[got[s][1]]] for s in sel])
            types.append((i, e)); pairs.append(pr)
print(json.dumps({'rho': rho, 's_min': round(smin, 4), 'types': len(types), 'of': 5 * k, 'sec': round(time.time() - t0)}), flush=True)
T = len(types); comp = np.zeros((T, T), bool); lim = (1 + eps) ** 2
for a in range(T):
    Pa = pairs[a]
    for b in range(a + 1, T):
        if types[a][0] == types[b][0]: continue
        Pb = pairs[b]; ok = np.ones((len(Pa), len(Pb)), bool)
        for s in range(2):
            for t in range(2):
                X = Pa[:, s]; Y = Pb[:, t]; D = (X * X).sum(1)[:, None] + (Y * Y).sum(1)[None] - 2 * X @ Y.T; ok &= D <= lim
        comp[a, b] = comp[b, a] = bool(ok.any())
print(json.dumps({'compatible_type_pairs': int(np.triu(comp, 1).sum()), 'type_pairs': int(T * (T - 1) / 2), 'sec': round(time.time() - t0)}), flush=True)
sol = Solver(name='cd153'); var = {t: n + 1 for n, t in enumerate(types)}
for i in range(k):
    cl = [var[t] for t in types if t[0] == i]
    if not cl: print(json.dumps({'rotation_unblockable': i})); sys.exit()
    sol.add_clause(cl)
for a in range(T):
    for b in range(a + 1, T):
        if types[a][0] != types[b][0] and not comp[a, b]: sol.add_clause([-var[types[a]], -var[types[b]]])
res = sol.solve()
print(json.dumps({'k': k, 'CSP': 'SAT (pairwise-compatible blocking choice exists)' if res else 'UNSAT (no X blocks all k, sampled)', 'sec': round(time.time() - t0)}), flush=True)
