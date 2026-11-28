# Two-fan idea on an outer body Omega: several good simplex fans with small critical zones (pairs in one cone at distance
# >= 1).  If the critical zones of fans i and j are incompatible (no unit pair from zone i and zone j with all 4 cross
# distances <= 1), every diameter-1 X ⊂ Omega is split by fan i or fan j.  Sampled test.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary, d2
from p6_cones import rot, SIMP
specarg, nfan, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
cthr = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0; kthr = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
nb = int(sys.argv[6]) if len(sys.argv) > 6 else 4000
spec = json.loads(open(specarg).read()) if specarg.endswith('.json') else json.loads(specarg)
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
S = np.concatenate([boundary(rho, P, cuts, nb, rng), interior(rho, P, cuts, nb // 5, rng)]); n = len(S)
D = np.empty((n, n), np.float32)
for i0 in range(0, n, 500): D[i0:i0 + 500] = d2(S[i0:i0 + 500], S)
def lab(z): return np.argmax((S - z[6:]) @ (SIMP @ rot(z[:6]).T).T, 1)
def phi(z):
    l = lab(z); m = 0.0
    for j in range(5):
        idx = np.flatnonzero(l == j)
        if len(idx) > 1: m = max(m, float(D[np.ix_(idx, idx)].max()))
    return np.sqrt(m)
fans = []
for k in range(nfan):
    z = np.r_[rng.uniform(-np.pi, np.pi, 6), np.zeros(4)]; f = phi(z); sig = np.r_[np.full(6, 0.3), np.full(4, 0.03)]
    for it in range(250):
        Z = z[None] + sig[None] * rng.normal(size=(8, 10)); fs = [phi(q) for q in Z]; o = int(np.argmin(fs))
        if fs[o] < f: z, f = Z[o], fs[o]
        else: sig *= 0.93
        if sig[0] < 1e-3: break
    l = lab(z); crit = []
    for j in range(5):
        idx = np.flatnonzero(l == j); sub = D[np.ix_(idx, idx)]; a, b = np.nonzero(np.triu(sub >= cthr * cthr, 1)); crit += list(zip(idx[a], idx[b]))
    fans.append((f, z, crit)); print(json.dumps({'fan': k, 'max_piece': round(f, 4), 'critical_pairs': len(crit)}), flush=True)
comp = np.zeros((nfan, nfan), bool)
for i in range(nfan):
    for j in range(i + 1, nfan):
        A = np.array(fans[i][2]) if fans[i][2] else np.zeros((0, 2), int); B = np.array(fans[j][2]) if fans[j][2] else np.zeros((0, 2), int)
        if len(A) == 0 or len(B) == 0: comp[i, j] = comp[j, i] = False; continue
        A = A[rng.choice(len(A), min(len(A), 1500), replace=False)]; B = B[rng.choice(len(B), min(len(B), 1500), replace=False)]
        ok = np.ones((len(A), len(B)), bool)
        for s in range(2):
            for t in range(2): ok &= D[np.ix_(A[:, s], B[:, t])] <= kthr * kthr
        comp[i, j] = comp[j, i] = bool(ok.any())
print(json.dumps({'fan_values': [round(f[0], 4) for f in fans], 'incompatible_fan_pairs': [(i, j) for i in range(nfan) for j in range(i + 1, nfan) if not comp[i, j]]}), flush=True)
# joint CSP: choose one critical pair per fan (fans with no critical pair always work -> certificate trivially)
from pysat.solvers import Solver
cp = [np.array(f[2]) if f[2] else np.zeros((0, 2), int) for f in fans]
if any(len(c) == 0 for c in cp): print(json.dumps({'joint': 'some fan has no sampled critical pair'})); sys.exit()
cp = [c[rng.choice(len(c), min(len(c), 600), replace=False)] for c in cp]
var = {}; nv = 0
for i, c in enumerate(cp):
    for p in range(len(c)): nv += 1; var[(i, p)] = nv
sol = Solver(name='cd153')
for i, c in enumerate(cp): sol.add_clause([var[(i, p)] for p in range(len(c))])
ncl = 0
for i in range(nfan):
    for j in range(i + 1, nfan):
        A, B = cp[i], cp[j]; ok = np.ones((len(A), len(B)), bool)
        for s_ in range(2):
            for t in range(2): ok &= D[np.ix_(A[:, s_], B[:, t])] <= kthr * kthr
        a, b = np.nonzero(~ok)
        for p, q in zip(a.tolist(), b.tolist()): sol.add_clause([-var[(i, p)], -var[(j, q)]]); ncl += 1
res = sol.solve()
print(json.dumps({'joint_CSP_fans': nfan, 'clauses': ncl, 'result': 'SAT (some X blocks all fans, sampled)' if res else 'UNSAT (every X split by one of the fans, sampled)'}), flush=True)
if res:
    m = sol.get_model(); chosen = [p for (i, p), v in var.items() if m[v - 1] > 0]
    print(json.dumps({'note': 'compatible blocking choice exists'}))
