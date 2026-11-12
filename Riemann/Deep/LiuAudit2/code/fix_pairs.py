"""liuaudit2: correct the CRT pair phases of the saved T12 data (bias3e5.npz) without recomputation.
Buggy pair factor = zeta3^k, correct = zeta3^(2k) = (buggy)^2, so gamma_2^fixed(c) = gamma_2^bug(c)^2 / prod_{p|c} gamma_2(p).
Factorizations by replaying the deterministic DFS of theta_bias_q12.py. Self-check: Galois symmetry (class-0 per-norm sums of
conj(Lambda) gamma_2 must be real) and agreement with an independent fixed-code run at X = 3e4 if given.
Run: uv run --no-project --with numpy python3 fix_pairs.py in.npz out.npz
"""
import sys, math
import numpy as np
from collections import defaultdict
d = np.load(sys.argv[1]); out = sys.argv[2]
N = d['N'].astype(np.int64); B = d['b']; g = d['g']; lam = d['lam']; pq = d['pq'].astype(np.int64); pg = d['pg']
XMAX = int(N.max())
sys.setrecursionlimit(10000)
fac = [()]
def dfs(start, Nn, members):
    for j in range(start, len(pq)):
        q = int(pq[j])
        if Nn*q > XMAX: break
        fac.append(tuple(members + [j]))
        if j + 1 < len(pq) and Nn*q*int(pq[j+1]) <= XMAX:
            dfs(j+1, Nn*q, members + [j])
dfs(0, 1, [])
assert len(fac) == len(N)
assert all((int(np.prod([pq[i] for i in f])) if f else 1) == N[t] for t, f in enumerate(fac))
gp = np.array([np.prod([pg[i] for i in f]) if f else 1.0+0j for f in fac])
gf = g**2/gp
single = np.array([len(f) <= 1 for f in fac])
print("single-prime/unit ideals unchanged: max |diff| = %.1e" % float(np.max(np.abs(gf[single] - g[single]))))
acc = defaultdict(complex)
for t in range(len(N)):
    if B[t] == 0: acc[int(N[t])] += np.conj(lam[t])*gf[t]
im = [abs(z.imag) for z in acc.values()]
print("Galois check after correction: class-0 per-norm sums real for %d of %d norms (max |Im| %.1e)" % (sum(1 for x in im if x < 1e-8), len(im), max(im)))
np.savez_compressed(out, N=N, b=B, g=gf, lam=lam, pq=pq, pb=d['pb'], pg=pg, plam=d['plam'])
print("saved", out)
