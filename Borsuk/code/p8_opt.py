# P8: rotated-cone optimisation against (nearly) true piece diameters: sampled max pair per piece refined by SLSQP
# local maximisation over the piece (balls, cuts, cone halfspaces); evolution strategy; then cert2 on all 5 cells.
import numpy as np, sys, json, subprocess
from scipy.optimize import minimize
sys.path.insert(0, 'Borsuk/code')
from p6_survey import random_config, boundary
from p6_cones import rot, SIMP, optimise
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
typ, seed, iters = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed); rho, P, cuts, lam = random_config(typ, rng)
X = boundary(rho, P, cuts, 8000, rng)
def piece_true(Q, ap, j, pairs=2):
    D = SIMP @ Q.T; lab = np.argmax(dot(X - ap, D), 1); Y = X[lab == j]
    if len(Y) < 2: return 0.0
    best = []
    for i0 in range(0, len(Y), 300):
        M = d2(Y[i0:i0 + 300], Y); k = np.argmax(M); a, b = np.unravel_index(k, M.shape); best.append((M[a, b], i0 + a, b))
    best.sort(reverse=True); val = best[0][0]
    cons = []
    for z in (0, 1):
        sl = slice(4 * z, 4 * z + 4)
        cons.append({'type': 'ineq', 'fun': lambda v, sl=sl: rho - v[sl] @ v[sl]})
        for p in P: cons.append({'type': 'ineq', 'fun': lambda v, sl=sl, p=p: 1 - ((v[sl] - p) ** 2).sum()})
        for n, c in cuts: cons.append({'type': 'ineq', 'fun': lambda v, sl=sl, n=n, c=c: c - v[sl] @ n})
        for k in range(5):
            if k != j: cons.append({'type': 'ineq', 'fun': lambda v, sl=sl, k=k: (v[sl] - ap) @ (D[j] - D[k])})
    for m, a, b in best[:pairs]:
        r = minimize(lambda v: -((v[:4] - v[4:]) ** 2).sum(), np.r_[Y[a], Y[b]], constraints=cons, method='SLSQP', options={'maxiter': 200, 'ftol': 1e-12})
        if r.success: val = max(val, -r.fun)
    return np.sqrt(val)
def obj(zv):
    Q = rot(zv[:6]); return max(piece_true(Q, zv[6:], j) for j in range(5))
f0, z = optimise(X[:3000], rng, iters=60)
mean = z.copy(); sig = np.array([0.05] * 6 + [0.01] * 4); best = (obj(z), z.copy())
for it in range(iters):
    Z = mean[None] + sig[None] * rng.normal(size=(8, 10)); f = np.array([obj(q) for q in Z]); o = np.argsort(f)
    if f[o[0]] < best[0]: best = (f[o[0]], Z[o[0]].copy())
    mean = 0.5 * mean + 0.5 * Z[o[:3]].mean(0); sig *= 0.97
    print(json.dumps({'it': it, 'best_true': round(best[0], 5)}), flush=True)
Q = rot(best[1][:6]); ap = best[1][6:]; D = SIMP @ Q.T; R = np.sqrt(rho)
txt = '5\n'
for j in range(5):
    hs = [(D[k] - D[j], (D[k] - D[j]) @ ap) for k in range(5) if k != j] + [(n, c) for n, c in cuts]
    txt += '%d %d\n0 0 0 0 %.18f\n' % (1 + len(P), len(hs), rho)
    for p in P: txt += ' '.join('%.18f' % x for x in p) + ' 1\n'
    for a, b in hs: txt += ' '.join('%.18f' % x for x in a) + ' %.18f\n' % b
txt += '5\n' + ''.join('%d %d 1.0 %.6f 3000000\n' % (j, j, R + 1e-9) for j in range(5))
open('/tmp/p8opt_%d_%d.txt' % (typ, seed), 'w').write(txt)
out = subprocess.run(['bash', '-c', 'Borsuk/.scratch/cert2 < /tmp/p8opt_%d_%d.txt' % (typ, seed)], capture_output=True, text=True).stdout
print(json.dumps({'type': typ, 'seed': seed, 'rho': round(rho, 4), 'best_true': round(best[0], 5), 'cert2': out.split('\n')[:5]}), flush=True)
