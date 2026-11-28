# P6/P8 relaxed test: 5-colour a dense boundary sample of an outer body so that no colour class has a pair at distance
# in the band [1 - eps, 1 + eps]; then check cluster consistency (each proximity-connected component of a colour
# class has diameter < 1).  If both hold on dense samples, partitions without unit pairs (unions of far clusters) exist.
import numpy as np, sys, json, time
sys.path.insert(0, 'Borsuk/code')
from tabu import colour_tabu
from p6_pd import body, interior
from p6_survey import boundary
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
spec = json.loads(sys.argv[1]); M = int(sys.argv[2]); eps = float(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
C = np.concatenate([boundary(rho, P, cuts, M, rng), interior(rho, P, cuts, M // 4, rng)]); M = len(C)
lo2, hi2 = (1 - eps) ** 2, (1 + eps) ** 2
U = []; V = []
for i0 in range(0, M, 100):
    D = d2(C[i0:i0 + 100], C); ii, jj = np.nonzero((D >= lo2) & (D <= hi2)); m = jj > ii + i0; U.append((ii[m] + i0).astype(np.int32)); V.append(jj[m].astype(np.int32))
U = np.concatenate(U); V = np.concatenate(V)
t0 = time.time(); best = None
for sd in range(1, 4):
    col, conf = colour_tabu(M, U, V, 5, iters=800000 * sd, seed=sd)
    if best is None or conf < best[1]: best = (col, conf)
    if conf == 0: break
col, conf = best
out = {'spec': spec, 'M': M, 'eps': eps, 'band edges': int(len(U)), 'tabu conflicts': conf, 'sec': round(time.time() - t0, 1)}
if conf == 0:
    # cluster consistency: proximity graph (distance < h) inside each colour class; component diameters
    h = 2.5 * (np.pi ** 2 * 2 * rho ** 1.5 / M) ** (1 / 3)
    worst = 0.0
    for c in range(5):
        idx = np.flatnonzero(col == c); Y = C[idx]; n = len(Y)
        parent = np.arange(n)
        def find(a):
            while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
            return a
        for i0 in range(0, n, 200):
            D = d2(Y[i0:i0 + 200], Y); ii, jj = np.nonzero(D < h * h)
            for a, b in zip((ii + i0).tolist(), jj.tolist()):
                ra, rb = find(a), find(b)
                if ra != rb: parent[ra] = rb
        roots = np.array([find(a) for a in range(n)])
        for r in np.unique(roots):
            Z = Y[roots == r]
            for i0 in range(0, len(Z), 200):
                if len(Z) > 1: worst = max(worst, float(d2(Z[i0:i0 + 200], Z).max()))
    out['link h'] = round(h, 3); out['max component diameter'] = round(np.sqrt(worst), 4)
print(json.dumps(out), flush=True)
