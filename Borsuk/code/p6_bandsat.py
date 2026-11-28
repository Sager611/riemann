import numpy as np, sys, json, subprocess
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
spec = json.loads(sys.argv[1]); M = int(sys.argv[2]); eps = float(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
C = np.concatenate([boundary(rho, P, cuts, M, rng), interior(rho, P, cuts, M // 4, rng)]); M = len(C)
lo2, hi2 = (1 - eps) ** 2, (1 + eps) ** 2; E = []
for i0 in range(0, M, 100):
    D = d2(C[i0:i0 + 100], C); ii, jj = np.nonzero((D >= lo2) & (D <= hi2)); m = jj > ii + i0; E.append(np.stack([ii[m] + i0, jj[m]], 1))
E = np.concatenate(E); k = 5
with open('/tmp/band2.cnf', 'w') as f:
    f.write('p cnf %d %d\n' % (M * k, M + k * len(E) + 1))
    for v in range(M): f.write(' '.join(str(v * k + c + 1) for c in range(k)) + ' 0\n')
    for a, b in E.tolist():
        for c in range(k): f.write('-%d -%d 0\n' % (a * k + c + 1, b * k + c + 1))
    f.write('%d 0\n' % (int(E[0, 0]) * k + 1))
r = subprocess.run(['Borsuk/.scratch/tools/kissat/build/kissat', '/tmp/band2.cnf', '/tmp/band2.drat', '--no-binary'], capture_output=True, text=True)
res = {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, r.returncode)
print(json.dumps({'spec': spec, 'M': M, 'eps': eps, 'edges': int(len(E)), 'kissat': res}), flush=True)
if res == 'UNSAT':
    v = subprocess.run(['Borsuk/.scratch/tools/drat-trim/drat-trim', '/tmp/band2.cnf', '/tmp/band2.drat'], capture_output=True, text=True).stdout
    print([l for l in v.splitlines() if l.startswith('s ')]); np.save('Borsuk/.scratch/band_core.npy', C)

if res == 'SAT':
    lits = [int(t) for l in r.stdout.splitlines() if l.startswith('v ') for t in l[2:].split()]
    pos = set(l for l in lits if l > 0)
    col = np.array([next(c for c in range(k) if v * k + c + 1 in pos) for v in range(M)])
    assert all(col[a] != col[b] for a, b in E.tolist())
    for h in (0.15, 0.2, 0.275):
        worst = 0.0
        for c in range(5):
            Y = C[col == c]; n = len(Y); parent = list(range(n))
            def find(a):
                while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
                return a
            for i0 in range(0, n, 200):
                D = d2(Y[i0:i0 + 200], Y); ii, jj = np.nonzero(D < h * h)
                for a, b in zip((ii + i0).tolist(), jj.tolist()):
                    ra, rb = find(a), find(b)
                    if ra != rb: parent[ra] = rb
            roots = np.array([find(a) for a in range(n)])
            for rt in np.unique(roots):
                Z = Y[roots == rt]
                if len(Z) > 1: worst = max(worst, max(float(d2(Z[i0:i0 + 200], Z).max()) for i0 in range(0, len(Z), 200)))
        print(json.dumps({'link h': h, 'max component diameter': round(np.sqrt(worst), 4)}), flush=True)
