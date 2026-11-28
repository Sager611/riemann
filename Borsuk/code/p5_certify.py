# Certify (non-)5-colourability of the far graph (pairs at distance >= thr) of a point set Y (float coords):
# 5-core + dominance reduction, clique-based symmetry breaking, kissat (+ DRAT check by drat-trim when UNSAT).
import numpy as np, sys, os, subprocess, json, time
thr = float(sys.argv[2]); Y = np.load(sys.argv[1]); n = len(Y); tag = sys.argv[3]
T = 'Borsuk/.scratch/tools'
def dot(A, B): return sum(A[:, k, None] * B[None, :, k] for k in range(A.shape[1]))
def d2(A, B): return (A * A).sum(1)[:, None] + (B * B).sum(1)[None] - 2 * dot(A, B)
Adj = np.zeros((n, n), bool)
for i0 in range(0, n, 100): Adj[i0:i0 + 100] = d2(Y[i0:i0 + 100], Y) >= thr * thr
np.fill_diagonal(Adj, False)
alive = np.ones(n, bool); changed = True
while changed:
    changed = False; idx = np.flatnonzero(alive); A = Adj[np.ix_(idx, idx)]; deg = A.sum(1)
    if (deg < 5).any(): alive[idx[deg < 5]] = False; changed = True; continue
    for a in range(len(idx)):
        if not alive[idx[a]]: continue
        Na = A[a] & alive[idx]; cand = (~A[a]) & alive[idx]; cand[a] = False
        if (~np.any(Na[None, :] & ~A[cand], axis=1)).any(): alive[idx[a]] = False; changed = True
idx = np.flatnonzero(alive); A = Adj[np.ix_(idx, idx)]; m = len(idx)
U, V = np.nonzero(np.triu(A, 1))
# greedy max clique for symmetry breaking
best = []
for s in np.argsort(-A.sum(1))[:200]:
    cl = [s]; cand = np.flatnonzero(A[s])
    while len(cand):
        v = cand[np.argmax(A[np.ix_(cand, cand)].sum(1))]; cl.append(v); cand = cand[A[v, cand]]
    if len(cl) > len(best): best = cl
k = 5; ncl = m + len(U) * k + len(best)
with open('/tmp/%s.cnf' % tag, 'w') as f:
    f.write('p cnf %d %d\n' % (m * k, m + len(U) * k + len(best)))
    for v in range(m): f.write(' '.join(str(v * k + c + 1) for c in range(k)) + ' 0\n')
    for a, b in zip(U.tolist(), V.tolist()):
        for c in range(k): f.write('-%d -%d 0\n' % (a * k + c + 1, b * k + c + 1))
    for c, v in enumerate(best[:k]): f.write('%d 0\n' % (v * k + c + 1))
print(json.dumps({'n': n, 'reduced': m, 'edges': int(len(U)), 'clique': len(best)}), flush=True)
t0 = time.time()
r = subprocess.run([T + '/kissat/build/kissat', '-q', '/tmp/%s.cnf' % tag, '/tmp/%s.drat' % tag, '--no-binary'], capture_output=True, text=True)
res = 'SAT' if r.returncode == 10 else 'UNSAT' if r.returncode == 20 else 'rc%d' % r.returncode
print(json.dumps({'kissat': res, 'sec': round(time.time() - t0, 1)}), flush=True)
if res == 'UNSAT':
    v = subprocess.run([T + '/drat-trim/drat-trim', '/tmp/%s.cnf' % tag, '/tmp/%s.drat' % tag], capture_output=True, text=True).stdout
    print('drat-trim:', [l for l in v.splitlines() if l.startswith('s ')])
    np.save('Borsuk/.scratch/%s_core.npy' % tag, Y[idx])
