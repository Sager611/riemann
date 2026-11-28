# Shrink a barrier point set: block deletion keeping the far graph (distance >= t) not 5-colourable (kissat, float
# distances with a safety margin: an edge is kept only if its float distance >= t + 1e-9).  Writes the reduced file.
import numpy as np, sys, subprocess, os, time
from fractions import Fraction as Fr
src, out, t = sys.argv[1], sys.argv[2], float(sys.argv[3])
lines = [l for l in open(src) if l.strip() and not l.startswith('#')]
X = np.array([[float(Fr(c)) for c in l.split()] for l in lines]); n = len(X)
D = np.empty((n, n), np.float32)
for i0 in range(0, n, 500):
    A = X[i0:i0 + 500]; D[i0:i0 + 500] = np.sqrt(np.maximum((A * A).sum(1)[:, None] + (X * X).sum(1)[None] - 2 * A @ X.T, 0))
E = D >= t + 1e-9
def unsat(keep):
    idx = np.flatnonzero(keep); sub = E[np.ix_(idx, idx)]; U, V = np.nonzero(np.triu(sub, 1)); m = len(idx); k = 5
    fn = '/tmp/shrink_%d.cnf' % os.getpid()
    with open(fn, 'w') as f:
        f.write('p cnf %d %d\n' % (m * k, m + k * len(U) + 2))
        np.savetxt(f, np.c_[(k * np.arange(m))[:, None] + np.arange(1, k + 1)[None, :], np.zeros(m, int)], fmt='%d')
        for c in range(k): np.savetxt(f, np.c_[-(k * U + c + 1), -(k * V + c + 1), np.zeros(len(U), int)], fmt='%d')
        f.write('%d 0\n%d 0\n' % (k * U[0] + 1, k * V[0] + 2))
    try: r = subprocess.run(['Borsuk/.scratch/tools/kissat/build/kissat', '-q', fn], capture_output=True, timeout=900).returncode
    except subprocess.TimeoutExpired: r = 0
    os.remove(fn); return r == 20
keep = np.ones(n, bool); rng = np.random.default_rng(0); t0 = time.time()
assert unsat(keep)
for blk in [300, 150, 80, 40, 20, 10]:
    idx = np.flatnonzero(keep); rng.shuffle(idx)
    for s in range(0, len(idx), blk):
        if time.time() - t0 > 10000: break
        trial = keep.copy(); trial[idx[s:s + blk]] = False
        if unsat(trial): keep = trial
    print('block', blk, 'kept', int(keep.sum()), 'sec', round(time.time() - t0), flush=True)
with open(out, 'w') as f:
    f.write(open(src).readline())
    for i in np.flatnonzero(keep): f.write(lines[i])
print('wrote', int(keep.sum()), 'points to', out)
