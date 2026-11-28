# KEY RELAXATION: a part of the outer body may contain pairs at distance > 1 (they never lie in one set of
# diameter 1); it must contain no pair at distance exactly 1.  With parts = unions of convex blobs: each blob has
# diameter < 1, and two blobs of one part are either 'close' (all cross distances < 1) or 'far' (all > 1), and
# close is transitive inside a part.  Sample-based test: blobs = Voronoi cells of K centres; SAT (kissat) for a
# 5-colouring with pair conflicts (interval [min-d, max+d] meets 1) and triple conflicts (A~B, B~C close, A,C far).
import numpy as np, sys, json, time, subprocess
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
spec = json.loads(sys.argv[1]); K = int(sys.argv[2]); delta = float(sys.argv[3]); seed = int(sys.argv[4]); NB = int(sys.argv[5])
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
S = np.concatenate([boundary(rho, P, cuts, NB, rng), interior(rho, P, cuts, NB // 3, rng)])
# k-means-ish centres (few Lloyd steps)
C = S[rng.choice(len(S), K, replace=False)]
for _ in range(8):
    lab = np.concatenate([np.argmin(d2(S[i:i + 2000], C), 1) for i in range(0, len(S), 2000)])
    C = np.array([S[lab == k].mean(0) if np.any(lab == k) else C[k] for k in range(K)])
lab = np.concatenate([np.argmin(d2(S[i:i + 2000], C), 1) for i in range(0, len(S), 2000)])
blobs = [S[lab == k] for k in range(K)]
mn = np.full((K, K), 9.0); mx = np.zeros((K, K))
for a in range(K):
    Da = [np.sqrt(np.maximum(d2(blobs[a], blobs[b]), 0)) if len(blobs[a]) and len(blobs[b]) else None for b in range(K)]
    for b in range(K):
        if Da[b] is not None and Da[b].size: mn[a, b] = Da[b].min(); mx[a, b] = Da[b].max()
diam = np.diag(mx).max(); close = mx + delta < 1; far = mn - delta > 1
np.fill_diagonal(close, True)
conf = ~(close | far)
print(json.dumps({'K': K, 'max blob diam': round(float(diam), 3), 'pair conflicts': int(np.triu(conf, 1).sum()), 'close pairs': int(np.triu(close, 1).sum()), 'far pairs': int(np.triu(far, 1).sum())}), flush=True)
k = 5; cl = []
for v in range(K): cl.append([v * k + c + 1 for c in range(k)])
I, J = np.nonzero(np.triu(conf, 1))
for a, b in zip(I.tolist(), J.tolist()):
    for c in range(k): cl.append([-(a * k + c + 1), -(b * k + c + 1)])
ntri = 0
for b in range(K):
    nb = np.flatnonzero(close[b]); nb = nb[nb != b]
    sub = far[np.ix_(nb, nb)]; ii, jj = np.nonzero(np.triu(sub, 1))
    for a, cc in zip(nb[ii].tolist(), nb[jj].tolist()):
        ntri += 1
        for c in range(k): cl.append([-(a * k + c + 1), -(b * k + c + 1), -(cc * k + c + 1)])
cl.append([1])
with open('/tmp/blob.cnf', 'w') as f:
    f.write('p cnf %d %d\n' % (K * k, len(cl)))
    for c_ in cl: f.write(' '.join(map(str, c_)) + ' 0\n')
t0 = time.time(); r = subprocess.run(['Borsuk/.scratch/tools/kissat/build/kissat', '-q', '/tmp/blob.cnf'], capture_output=True, text=True)
print(json.dumps({'triples': ntri, 'clauses': len(cl), 'kissat': {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, r.returncode), 'sec': round(time.time() - t0, 1)}), flush=True)
