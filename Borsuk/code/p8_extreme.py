# Extreme all-contact bodies: type 3 equilateral triangle at rho = 1/3 (slab cuts), type 4 regular tetrahedron at rho = 3/8
# (cut x4 <= 1/2), slightly smaller rho variants (scaled triangle/tetrahedron). CEGAR beta5 bracket on boundary samples.
import numpy as np, sys, json
sys.path.insert(0, 'Borsuk/code')
from p6_survey import boundary, cegar
which, rho, M, seed = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
thrs = [float(x) for x in sys.argv[5].split(',')]
R = np.sqrt(rho); rng = np.random.default_rng(seed)
if which == '3':
    th = 2 * np.pi * np.arange(3) / 3; P = R * np.stack([np.cos(th), np.sin(th), 0 * th, 0 * th], 1)
    cuts = [(np.array([0, 0, 1.0, 0]), 0.5), (np.array([0, 0, -1.0, 0]), 0.5), (np.array([0, 0, 0, 1.0]), 0.5)]
else:
    v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3); P = R * np.hstack([v, np.zeros((4, 1))])
    cuts = [(np.array([0, 0, 0, 1.0]), 0.5)]
C = boundary(rho, P, cuts, M, rng); res = {}
for t in thrs:
    r = cegar(C, t, TL=900, rng=rng); res[t] = r
    print(json.dumps({'type': which, 'rho': rho, 'maxpair': round(float(np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1).max())), 4), 'thr': t, 'res': r}), flush=True)
    if r[0] == 'COL': break
