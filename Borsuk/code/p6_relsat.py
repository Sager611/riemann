# Cluster-consistent relaxed colouring: colour c is single-cluster (no same-colour pair >= 1-eps) or 'double' on a
# far region pair (A_j, B_j) (all cross distances > 1+eps): then its points lie in A_j ∪ B_j and only A_j × B_j pairs
# are exempt.  kissat.
import numpy as np, sys, json, subprocess, itertools
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
spec = json.loads(sys.argv[1]); M = int(sys.argv[2]); eps = float(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
if len(sys.argv) > 5:
    from fractions import Fraction as Fr
    C = np.array([[float(Fr(t)) for t in l.split()] for l in open(sys.argv[5]) if l.strip() and not l.startswith('#')])
else:
    C = np.concatenate([boundary(rho, P, cuts, M, rng), interior(rho, P, cuts, M // 4, rng)])
M = len(C)
D = np.empty((M, M), np.float32)
for i0 in range(0, M, 200): D[i0:i0 + 200] = np.sqrt(np.maximum(d2(C[i0:i0 + 200], C), 0))
# candidate far region pairs from directions
dirs = []
for v in itertools.product((-1, 0, 1), repeat=4):
    v = np.array(v, float)
    if np.any(v) and tuple(-v) not in [tuple(x) for x in dirs]: dirs.append(v / np.linalg.norm(v))
pairs = []
for w in dirs:
    h = C @ w
    for t in np.linspace(0.30, 0.60, 13):
        A = np.flatnonzero(h >= t); B = np.flatnonzero(h <= -t)
        if len(A) < 5 or len(B) < 5: continue
        if D[np.ix_(A, B)].min() > 1 + eps: pairs.append((A, B)); break
k = 5; J = len(pairs); nv = M * k
var_D = lambda c, j: nv + c * J + j + 1
cl = []
for v in range(M): cl.append([v * k + c + 1 for c in range(k)])
exempt = np.zeros((M, M), bool)
inreg = np.zeros((J, M), bool)
for j, (A, B) in enumerate(pairs):
    exempt[np.ix_(A, B)] = True; exempt[np.ix_(B, A)] = True; inreg[j, A] = True; inreg[j, B] = True
I, Jj = np.nonzero(np.triu(D >= 1 - eps, 1))
nex = 0
for a, b in zip(I.tolist(), Jj.tolist()):
    if exempt[a, b]:
        # same colour allowed only if that colour is double on a region pair containing a in A_j and b in B_j
        js = [j for j in range(J) if inreg[j, a] and inreg[j, b] and ((a in set(pairs[j][0]) and b in set(pairs[j][1])) or (a in set(pairs[j][1]) and b in set(pairs[j][0])))]
        for c in range(k): cl.append([-(a * k + c + 1), -(b * k + c + 1)] + [var_D(c, j) for j in js])
        nex += 1
    else:
        for c in range(k): cl.append([-(a * k + c + 1), -(b * k + c + 1)])
for c in range(k):
    for j in range(J):
        for v in np.flatnonzero(~inreg[j]).tolist(): cl.append([-var_D(c, j), -(v * k + c + 1)])
    for j1 in range(J):
        for j2 in range(j1 + 1, J): cl.append([-var_D(c, j1), -var_D(c, j2)])
with open('/tmp/rel.cnf', 'w') as f:
    f.write('p cnf %d %d\n' % (nv + k * J, len(cl)))
    for c_ in cl: f.write(' '.join(map(str, c_)) + ' 0\n')
r = subprocess.run(['Borsuk/.scratch/tools/kissat/build/kissat', '-q', '/tmp/rel.cnf'], capture_output=True, text=True)
print(json.dumps({'spec': spec, 'M': M, 'eps': eps, 'region pairs': J, 'exempt far pairs': nex, 'clauses': len(cl), 'kissat': {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, r.returncode)}), flush=True)
