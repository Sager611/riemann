# Cluster-consistent relaxed 5-colouring (as p6_relsat.py) with streamed CNF: colour c is single-cluster (no same-colour
# pair >= 1-eps) or 'double' on one far region pair (A_j, B_j) (all cross distances > 1+eps; colour-c points inside A_j ∪ B_j,
# only A_j x B_j pairs exempt).  Also runs the plain encoding on the same sample for comparison.
import numpy as np, sys, json, subprocess, itertools, os
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary, d2
spec = json.loads(open(sys.argv[1]).read()) if sys.argv[1].endswith('.json') else json.loads(sys.argv[1])
M = int(sys.argv[2]); eps = float(sys.argv[3]); seed = int(sys.argv[4]); mode = sys.argv[5] if len(sys.argv) > 5 else 'both'
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
C = np.concatenate([boundary(rho, P, cuts, M, rng), interior(rho, P, cuts, M // 4, rng)]); N = len(C)
lo2, hi2 = (1 - eps) ** 2, (1 + eps) ** 2
I, J = [], []
for i0 in range(0, N, 500):
    D = d2(C[i0:i0 + 500], C); a, b = np.nonzero(D >= lo2); m = b > a + i0; I.append(a[m] + i0); J.append(b[m])
I = np.concatenate(I).astype(np.int64); J = np.concatenate(J).astype(np.int64)
dirs = []
for v in itertools.product((-1, 0, 1), repeat=4):
    v = np.array(v, float)
    if np.any(v) and not any(np.allclose(-v / np.linalg.norm(v), w) for w in dirs): dirs.append(v / np.linalg.norm(v))
regs = []
for w in dirs:
    h = C @ w
    for t in np.linspace(0.30, 0.60, 13):
        A = h >= t; B = h <= -t
        if A.sum() < 5 or B.sum() < 5: continue
        ok = True
        for i0 in range(0, N, 500):
            Ai = np.flatnonzero(A[i0:i0 + 500]) + i0
            if len(Ai) and d2(C[Ai], C[B]).min() <= hi2: ok = False; break
        if ok: regs.append((A, B)); break
R = len(regs); assert R <= 62
mA = np.zeros(N, np.int64); mB = np.zeros(N, np.int64)
for j, (A, B) in enumerate(regs): mA[A] |= 1 << j; mB[B] |= 1 << j
def run(relaxed):
    k = 5; nv = N * k; vD = lambda c, j: nv + c * R + j + 1
    fn = '/tmp/rel2_%d.cnf' % os.getpid()
    js = ((mA[I] & mB[J]) | (mB[I] & mA[J])) if relaxed else np.zeros(len(I), np.int64)
    lines = 0
    with open(fn, 'w') as f:
        f.write('p cnf 0 0\n')  # placeholder, rewritten below
        base = (k * np.arange(N))[:, None] + np.arange(1, k + 1)[None, :]
        np.savetxt(f, np.c_[base, np.zeros(N, np.int64)], fmt='%d'); lines += N
        plain = js == 0
        for c in range(k):
            np.savetxt(f, np.c_[-(k * I[plain] + c + 1), -(k * J[plain] + c + 1), np.zeros(plain.sum(), np.int64)], fmt='%d'); lines += int(plain.sum())
        if relaxed:
            for val in np.unique(js[~plain]):
                sel = np.flatnonzero(js == val); bits = [j for j in range(R) if (int(val) >> j) & 1]
                for c in range(k):
                    extra = np.tile(np.array([vD(c, j) for j in bits], np.int64), (len(sel), 1))
                    np.savetxt(f, np.c_[-(k * I[sel] + c + 1), -(k * J[sel] + c + 1), extra, np.zeros(len(sel), np.int64)], fmt='%d'); lines += len(sel)
            for c in range(k):
                for j, (A, B) in enumerate(regs):
                    out = np.flatnonzero(~(A | B))
                    np.savetxt(f, np.c_[np.full(len(out), -vD(c, j)), -(k * out + c + 1), np.zeros(len(out), np.int64)], fmt='%d'); lines += len(out)
                for j1 in range(R):
                    for j2 in range(j1 + 1, R): f.write('%d %d 0\n' % (-vD(c, j1), -vD(c, j2))); lines += 1
    with open(fn, 'r+') as f:
        head = f.readline(); rest = f.read()
    with open(fn, 'w') as f:
        f.write('p cnf %d %d\n' % (nv + k * R, lines)); f.write(rest)
    r = subprocess.run(['Borsuk/.scratch/tools/kissat/build/kissat', '-q', fn], capture_output=True, text=True, timeout=7200)
    os.remove(fn)
    return {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, r.returncode)
out = {'N': N, 'eps': eps, 'pairs': int(len(I)), 'regions': R, 'exempt': int(((mA[I] & mB[J]) | (mB[I] & mA[J]) != 0).sum())}
if mode in ('both', 'plain'): out['plain'] = run(False)
if mode in ('both', 'relaxed'): out['relaxed'] = run(True)
print(json.dumps(out), flush=True)
