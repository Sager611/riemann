# torus_sos2.py: vectorized version of torus_sos.py (sparse linear map on vec(G)).
import sys, math, itertools, numpy as np, scipy.sparse as sp, cvxpy as cp, time
c = 10.0; L = math.log(c)
fac = float(sys.argv[1]); mode = sys.argv[2]; D2 = int(sys.argv[3]); D5 = int(sys.argv[4])
P = 2 * L * fac
pr = [2, 3, 5, 7]; lg = np.log(pr)
box = [range(-D2, D2 + 1), range(-D2, D2 + 1), range(-D5, D5 + 1), range(-D5, D5 + 1)]
B = np.array(list(itertools.product(*box))); nB = len(B)
diffs = (B[:, None, :] - B[None, :, :]).reshape(-1, 4)
keys, inv = np.unique(diffs, axis=0, return_inverse=True)
inv = inv.reshape(-1)
nK = len(keys)
A = sp.csr_matrix((np.ones(nB * nB), (inv, np.arange(nB * nB))), shape=(nK, nB * nB))
rows_eq = []; rhs = []; one_row = None
for r_, d in enumerate(keys):
    nz = np.nonzero(d)[0]
    lv = abs(float(np.dot(d, lg)))
    if len(nz) == 0: one_row = r_; continue
    if len(nz) == 1 and lv < L - 1e-12:
        p = pr[nz[0]]; k = abs(int(d[nz[0]]))
        rows_eq.append(r_); rhs.append(-math.log(p) / p**(k / 2)); continue
    gapok = (L - 1e-12 <= (lv % P) <= P - L + 1e-12)
    if mode == 'sep' and len(nz) > 1: gapok = False
    if mode == 'none': gapok = False
    if not gapok: rows_eq.append(r_); rhs.append(0.0)
Aeq = A[rows_eq, :]; a1 = A[one_row, :]
G = cp.Variable((nB, nB), symmetric=True)
beta = cp.Variable()
g = cp.vec(G, order='C')
prob = cp.Problem(cp.Minimize(beta), [G >> 0, Aeq @ g == np.array(rhs), a1 @ g == beta])
t0 = time.time()
prob.solve(solver=cp.CLARABEL)
sig2 = 2 * sum(math.log(p) / p**(k / 2) for p in pr for k in range(1, 10) if p**k < c)
nfree = nK - 1 - len(rows_eq)
print(f'c=10 P={fac}*2L mode={mode} box=({D2},{D2},{D5},{D5}): |B|={nB}, free atoms={nfree}, {prob.status}, envelope 2*beta = {2*beta.value:.4f} (2sigma {sig2:.4f}; ||C_p||_window 2.510) [{time.time()-t0:.1f}s]', flush=True)
