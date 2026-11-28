# Where is the 5-colouring obstruction of the extreme type-3 body?  CEGAR at thr; reduce the NOT set to a 5-core
# (iteratively drop points of degree < 5); kissat UNSAT check of the core; statistics of |x|, |a|, |y|, which boundary part.
import numpy as np, sys, json, subprocess, os
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, cegar, d2
rho, M, seed, thr = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
R = np.sqrt(rho); rng = np.random.default_rng(seed)
th = 2 * np.pi * np.arange(3) / 3; P = R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1)
cuts = [(np.array([0, 0, 1.0, 0]), 0.5), (np.array([0, 0, -1.0, 0]), 0.5), (np.array([0, 0, 0, 1.0]), 0.5)]
C = boundary(rho, P, cuts, M, rng)
r = cegar(C, thr, TL=900, rng=rng, want=True)
print(r[0], r[1], flush=True)
if r[0] != 'NOT': sys.exit()
X = r[2]; N = len(X)
E = d2(X, X) >= thr * thr; np.fill_diagonal(E, False)
alive = np.ones(N, bool)
while True:
    deg = (E & alive[None, :]).sum(1); drop = alive & (deg < 5)
    if not drop.any(): break
    alive &= ~drop
Y = X[alive]; n = len(Y); Ey = E[np.ix_(alive, alive)]
U, V = np.nonzero(np.triu(Ey, 1))
cnf = []
for i in range(n): cnf.append(' '.join(str(5 * i + c + 1) for c in range(5)) + ' 0')
for u, v in zip(U, V):
    for c in range(5): cnf.append('-%d -%d 0' % (5 * u + c + 1, 5 * v + c + 1))
# symmetry breaking: first edge's endpoints colours 1, 2
cnf.append('%d 0' % (5 * U[0] + 1)); cnf.append('%d 0' % (5 * V[0] + 2))
cnf.insert(0, 'p cnf %d %d' % (5 * n, len(cnf)))
open('/tmp/ext3core.cnf', 'w').write('\n'.join(cnf) + '\n')
k = subprocess.run(['Borsuk/.scratch/tools/kissat/build/kissat', '-q', '/tmp/ext3core.cnf'], capture_output=True, text=True, timeout=3000)
nx = np.linalg.norm(Y, axis=1); na = np.linalg.norm(Y[:, :2], axis=1); ny = np.linalg.norm(Y[:, 2:], axis=1)
part = np.where(nx > R - 1e-9, 'sphere', 'other')
dk = np.sqrt(((Y[:, None, :] - P[None]) ** 2).sum(-1)).max(1)
part = np.where(np.abs(dk - 1) < 1e-9, 'ball_p', part)
part = np.where((np.abs(np.abs(Y[:, 2]) - 0.5) < 1e-9) | (np.abs(Y[:, 3] - 0.5) < 1e-9), 'cut', part)
vals, cnts = np.unique(part, return_counts=True)
print(json.dumps({'core': n, 'edges': len(U), 'kissat_rc': k.returncode, 'parts': dict(zip(vals.tolist(), cnts.tolist())),
                  'norm_hist': np.histogram(nx, bins=8, range=(0.3, R + 1e-6))[0].tolist(),
                  'na_hist': np.histogram(na, bins=8, range=(0, R + 1e-6))[0].tolist(), 'ny_hist': np.histogram(ny, bins=8, range=(0, R + 1e-6))[0].tolist()}), flush=True)
np.save('Borsuk/.scratch/ext3_core.npy', Y)
