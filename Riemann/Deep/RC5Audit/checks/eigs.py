# eigs.py (rc5audit): top eigenvalues of the family Gram operator on admissible n in [N, 2N).
# usage: eigs.py prefix X N [Q0]
import sys, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh
pre, X, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); Q0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
cls = np.loadtxt(pre + '.cls', dtype=np.int64, ndmin=2)  # D h ng
trip = np.fromfile(pre + '.trip', dtype=np.int32).reshape(-1, 3)
gen = np.fromfile(pre + '.gen', dtype=np.int8)
def admissible(lo, hi):
    a = np.ones(hi, bool); a[0] = False
    for p in range(2, int(hi ** 0.5) + 2):
        a[p * p::p * p] = False
    for p in (2, 3, 5, 7, 11, 13, 17):
        a[::p] = False
    a[:lo] = False
    return np.nonzero(a)[0]
ns = admissible(N, 2 * N); m = len(ns); idx = -np.ones(2 * N, np.int64); idx[ns] = np.arange(m)
nc = int(cls[:, 1].sum())
# trip rows are written in class order: build CSR directly (low memory)
assert np.all(np.diff(trip[:, 0]) >= 0)
indptr = np.searchsorted(trip[:, 0], np.arange(nc + 1)).astype(np.int64)
indices = idx[trip[:, 1]].astype(np.int32); data = trip[:, 2].astype(np.float64) / 2.0
rowsV = trip[:, 0].copy(); del trip
V = sp.csr_matrix((data, indices, indptr), shape=(nc, m)); V.sum_duplicates()
hcl = np.repeat(cls[:, 1].astype(float), cls[:, 1])
# genus matrix: rows (D, g), cols classes
rows, cols, vals = [], [], []; r0 = 0; c0 = 0; gpos = 0
for D, h, ng in cls:
    blk = gen[gpos:gpos + h * ng].reshape(h, ng).astype(float); gpos += h * ng
    for g in range(ng):
        rows.append(np.full(h, r0 + g)); cols.append(np.arange(c0, c0 + h)); vals.append(blk[:, g])
    r0 += ng; c0 += h
G = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(r0, nc))
T = sp.csr_matrix((np.ones(nc), (np.repeat(np.arange(len(cls)), cls[:, 1]), np.arange(nc))), shape=(len(cls), nc))
# validation: sum over classes of D of v(n) = prod_{p|n}(1 + (D/p))
def kron(D, p):
    r = pow(D % p, (p - 1) // 2, p); return 0 if r == 0 else (1 if r == 1 else -1)
def pf(n):
    out = []; d = 19
    while d * d <= n:
        if n % d == 0: out.append(d); n //= d
        d += 2
    if n > 1: out.append(n)
    return out
TV = (T @ V).tocsc(); bad = 0; rng = np.random.default_rng(1)
for j in rng.choice(m, min(m, 40), replace=False):
    n = int(ns[j]); col = TV[:, j].toarray().ravel()
    for i in rng.choice(len(cls), min(len(cls), 50), replace=False):
        D = int(cls[i, 0]); exp = 1
        for p in pf(n): exp *= 1 + kron(D, p)
        if abs(col[i] - exp) > 1e-9: bad += 1
fam = hcl.sum() / 1.0  # = sum_D h(D)^2? no: hcl repeats h per class -> sum_D h^2; family size below
famsize = cls[:, 1].sum()
def op(kind, P=None):
    def mv(b):
        b = np.asarray(b).ravel()
        if P is not None: b = b - P @ (P.T @ b)
        y = V @ b; z = hcl * y
        if kind == 'notriv': z = z - T.T @ (T @ y)
        if kind == 'nongenus': z = z - G.T @ (G @ y)
        out = V.T @ z
        if P is not None: out = out - P @ (P.T @ out)
        return out
    return LinearOperator((m, m), matvec=mv, dtype=float)
res = {}
for kind in ('full', 'notriv', 'nongenus'):
    w = eigsh(op(kind), k=3, which='LA', tol=1e-5, return_eigenvectors=False)
    res[kind] = sorted(w, reverse=True)
Vc = V.tocoo(); diag = np.bincount(Vc.col, weights=Vc.data ** 2 * hcl[Vc.row], minlength=m); del Vc
line = f"X={X} N={N} m={m} famsize={famsize} validation_bad={bad} meandiag/fam={diag.mean()/famsize:.3f} maxdiag/fam={diag.max()/famsize:.3f} nongenus/maxdiag={res['nongenus'][0]/diag.max():.3f}"
for kind in res:
    line += f" | {kind}: " + " ".join(f"{x/famsize:.3f}" for x in res[kind])
print(line)
if Q0:
    # project out u_{q,psi}(n) = psi(n)(n/q), q odd squarefree <= Q0, psi mod 8
    vecs = []
    def jac(a, n):  # Jacobi (a/n), n odd positive
        a %= n; r = 1
        while a:
            while a % 2 == 0:
                a //= 2
                if n % 8 in (3, 5): r = -r
            a, n = n, a
            if a % 4 == 3 and n % 4 == 3: r = -r
            a %= n
        return r if n == 1 else 0
    psis = [lambda n: 1, lambda n: 1 if n % 4 == 1 else -1, lambda n: 1 if n % 8 in (1, 7) else -1, lambda n: 1 if n % 8 in (1, 3) else -1]
    for q in range(1, Q0 + 1, 2):
        if any(q % (p * p) == 0 for p in range(3, int(q ** 0.5) + 1, 2)): continue
        chiq = np.array([jac(int(n), q) for n in ns], float) if q > 1 else np.ones(m)
        for ps in psis:
            vecs.append(chiq * np.array([ps(int(n)) for n in ns], float))
    U = np.array(vecs).T; Qm, _ = np.linalg.qr(U)
    for kind in ('full', 'nongenus'):
        w = eigsh(op(kind, Qm), k=3, which='LA', tol=1e-5, return_eigenvectors=False)
        print(f"  projected out {U.shape[1]} quad-char dirs (q<={Q0}): {kind}: " + " ".join(f"{x/famsize:.3f}" for x in sorted(w, reverse=True)))
