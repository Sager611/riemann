# P11 feasibility probe (sampled, not rigorous): branching tree for the extreme type-3 body with X-dependent simplex fans.
# Node = sample of the leaf body.  If the best fan gives max piece diameter < 1 - margin: leaf.  Else take the worst pair
# (x, y) and branch on the cell C = B(x, delta): 'avoid' removes C; 'meet' intersects with B(x, 1 + delta).
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary, d2
from p6_cones import rot, SIMP
spec = json.loads(open(sys.argv[1]).read()); M = int(sys.argv[2]); delta = float(sys.argv[3]); margin = float(sys.argv[4]); budget = int(sys.argv[5])
rng = np.random.default_rng(1); rho, P, cuts = body(spec)
S0 = np.concatenate([boundary(rho, P, cuts, M, rng), interior(rho, P, cuts, M // 3, rng)])
def fan(S, z0=None):
    def obj(z):
        lab = np.argmax((S - z[6:]) @ (SIMP @ rot(z[:6]).T).T, 1); m = 0.0; arg = None
        for j in range(5):
            Y = S[lab == j]
            if len(Y) < 2: continue
            for i0 in range(0, len(Y), 1500):
                D = d2(Y[i0:i0 + 1500], Y); v = float(D.max())
                if v > m: m = v; a, b = np.unravel_index(int(np.argmax(D)), D.shape); arg = (Y[i0 + a], Y[b])
        return np.sqrt(m), arg
    best = (9, None, None)
    starts = [z0] if z0 is not None else []
    starts += [np.r_[rng.uniform(-np.pi, np.pi, 6), np.zeros(4)] for _ in range(2)]
    for z in starts:
        f, arg = obj(z); sig = np.r_[np.full(6, 0.2), np.full(4, 0.02)]
        for it in range(120):
            Z = z[None] + sig[None] * rng.normal(size=(6, 10)); out = [obj(q) for q in Z]; o = int(np.argmin([u[0] for u in out]))
            if out[o][0] < f: z, (f, arg) = Z[o], out[o]
            else: sig *= 0.9
            if sig[0] < 2e-3: break
        if f < best[0]: best = (f, z, arg)
    return best
t0 = time.time(); stack = [(np.ones(len(S0), bool), 0, None)]; nodes = leaves = 0; maxdepth = 0; worst_leaf = 0
while stack and nodes < budget:
    keep, depth, z0 = stack.pop(); nodes += 1; maxdepth = max(maxdepth, depth); S = S0[keep]
    f, z, arg = fan(S, z0)
    if f < 1 - margin or arg is None:
        leaves += 1; worst_leaf = max(worst_leaf, f)
    else:
        x = arg[0]; inC = ((S0 - x) ** 2).sum(1) <= delta ** 2; near = ((S0 - x) ** 2).sum(1) <= (1 + delta) ** 2
        stack.append((keep & ~inC, depth + 1, z)); stack.append((keep & near, depth + 1, z))
    if nodes % 20 == 0: print(json.dumps({'nodes': nodes, 'leaves': leaves, 'open': len(stack), 'maxdepth': maxdepth, 'last_f': round(float(f), 4), 'sec': round(time.time() - t0)}), flush=True)
print(json.dumps({'done': not stack, 'nodes': nodes, 'leaves': leaves, 'open': len(stack), 'maxdepth': maxdepth, 'worst_leaf': round(float(worst_leaf), 4), 'sec': round(time.time() - t0)}), flush=True)
