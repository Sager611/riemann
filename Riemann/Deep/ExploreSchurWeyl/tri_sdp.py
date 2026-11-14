#!/usr/bin/env python3
"""3-point SDP relaxation (explore50): tri_sym.py's LP plus the pivot positivity
  E sum_j (c + sum_k u(x_j - x_k))^2 >= 0 for all c, u  <=>  [[N, nu2^T], [nu2, nu3]] PSD  (nu3 as an M x M matrix).
Solved with cvxpy/Clarabel.  usage: tri_sdp.py N M [tri|pair] [CLARABEL|SCS]"""
import sys, numpy as np, cvxpy as cp
from scipy.sparse import csr_matrix
N = int(sys.argv[1]); M = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'tri'
gen = mode.endswith('g'); mode = mode.rstrip('g')     # 'trig'/'pairg': any multiplicities (m(m-1)(m-2) >= 0, no upper line bound)
X, Y = np.meshgrid(np.arange(M), np.arange(M), indexing='ij'); X = X.ravel(); Y = Y.ravel()
ims = []
for (u, v) in [(X, Y), (Y, X), (-X, Y-X), (Y-X, -X), (X-Y, -Y), (-Y, X-Y)]:
    for s in (1, -1): ims.append(((s*u) % M)*M + ((s*v) % M))
key = np.min(np.stack(ims), axis=0); uk, oid = np.unique(key, return_inverse=True); no = len(uk)
Bm = csr_matrix((np.ones(M*M), (np.arange(M*M), oid)), shape=(M*M, no))       # vec(nu3) = Bm v
ct = np.cos(2*np.pi*np.arange(M)/M); fr = np.arange(-(M//2) + 1, M//2 + 1)
rows, tg = [], []
for a in fr:
    for b in fr:
        if (abs(a) + abs(b) + abs(a+b))/2 > N: continue
        if mode == 'pair' and not (a == 0 and b == 0): continue
        t = N**3 if (a == 0 and b == 0) else (N*max(abs(a), abs(b)) if (a == 0 or b == 0 or a + b == 0) else 0.0)
        rows.append(np.bincount(oid, weights=ct[(a*X + b*Y) % M], minlength=no)); tg.append(t)
R = np.round(np.array(rows), 9); _, ui = np.unique(np.column_stack([R, np.array(tg)]), axis=0, return_index=True)
A3 = np.array(rows)[np.sort(ui)]; t3 = np.array(tg)[np.sort(ui)]
ks = np.arange(0, M//2 + 1); C2 = ct[np.outer(ks, np.arange(M)) % M]
v = cp.Variable(no, nonneg=True); n2 = cp.Variable(M, nonneg=True)
V = cp.reshape(Bm @ v, (M, M), order='C')
cons = [A3 @ v == t3, C2[:N+1] @ n2 == np.concatenate([[N**2], np.arange(1, N+1)]), C2[N+1:] @ n2 >= 0,
        cp.sum(V, axis=1) == N*n2]
line = V[0, :]
if gen: cons += [line >= n2, cp.sum(line) == N*n2[0], V[0, 0] >= 3*n2[0] - 2*N]
else: cons += [line >= n2, line <= 2*n2, cp.sum(line) == N*n2[0], V[0, 0] == 3*n2[0] - 2*N]
Z = cp.bmat([[np.array([[float(N)]]), cp.reshape(n2, (1, M), order='C')], [cp.reshape(n2, (M, 1), order='C'), V]])
cons += [(Z + Z.T)/2 >> 0]
pr = cp.Problem(cp.Maximize(n2[0]/N), cons)
slv = sys.argv[4] if len(sys.argv) > 4 else 'CLARABEL'
if slv == 'SCS': pr.solve(solver=cp.SCS, eps=1e-7, max_iters=200000)
else: pr.solve(solver=cp.CLARABEL)
print(f"N={N} M={M} mode={mode}{'(any m)' if gen else ''} SDP orbits={no} status={pr.status} maxM2={pr.value:.6f}", flush=True)
