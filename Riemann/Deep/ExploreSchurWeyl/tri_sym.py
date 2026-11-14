#!/usr/bin/env python3
"""Symmetry-reduced 3-point relaxation (explore50), same LP as tri_relax.py but with nu3 constant on orbits of the
order-12 group (S3 on the triple (j,k,l) x reflection) acting on Z_M^2; dedupes Fourier rows. Scales to N ~ 20, M ~ 12N.
usage: tri_sym.py N M [tri|pair]"""
import sys, numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
N = int(sys.argv[1]); M = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'tri'
X, Y = np.meshgrid(np.arange(M), np.arange(M), indexing='ij'); X = X.ravel(); Y = Y.ravel()
ims = []
for (u, v) in [(X, Y), (Y, X), (-X, Y-X), (Y-X, -X), (X-Y, -Y), (-Y, X-Y)]:
    for s in (1, -1):
        ims.append(((s*u) % M)*M + ((s*v) % M))
key = np.min(np.stack(ims), axis=0)
uk, oid = np.unique(key, return_inverse=True); no = len(uk)
ct = np.cos(2*np.pi*np.arange(M)/M)
fr = np.arange(-(M//2) + 1, M//2 + 1)
rows, tg = [], []
for a in fr:
    for b in fr:
        if (abs(a) + abs(b) + abs(a+b))/2 > N: continue
        if mode == 'pair' and not (a == 0 and b == 0): continue
        if a == 0 and b == 0: t = N**3
        elif a == 0 or b == 0 or a + b == 0: t = N*max(abs(a), abs(b))
        else: t = 0.0
        rows.append(np.bincount(oid, weights=ct[(a*X + b*Y) % M], minlength=no)); tg.append(t)
R = np.round(np.array(rows), 9); _, ui = np.unique(np.column_stack([R, np.array(tg)]), axis=0, return_index=True)
A3 = np.array(rows)[np.sort(ui)]; t3 = np.array(tg)[np.sort(ui)]
nv = M + no; off = M                                     # vars: nu2(0..M-1), v_orbit
ks = np.arange(0, M//2 + 1); C2 = ct[np.outer(ks, np.arange(M)) % M]
Aeq = [np.hstack([np.zeros((len(t3), M)), A3]), np.hstack([C2[:N+1], np.zeros((N+1, no))])]
beq = [t3, np.concatenate([[N**2], np.arange(1, N+1)])]
Aub = [np.hstack([-C2[N+1:], np.zeros((len(ks)-N-1, no))])]; bub = [np.zeros(len(ks)-N-1)]
Mg = np.zeros((M, nv))                                   # sum_y nu3(x,y) = N nu2(x)
for x in range(M):
    Mg[x, x] = -N; Mg[x, off:] = np.bincount(oid[x*M:(x+1)*M], minlength=no)
Aeq.append(Mg); beq.append(np.zeros(M))
lo = np.zeros((M, nv)); hi = np.zeros((M, nv)); sm = np.zeros((1, nv))
for t in range(M):                                       # line x = 0: nu3(0,t) = v_{o(0,t)}
    o = oid[0*M + t]
    lo[t, t] += 1.0; lo[t, off+o] -= 1.0; hi[t, off+o] += 1.0; hi[t, t] -= 2.0; sm[0, off+o] += 1.0
sm[0, 0] -= N
Aub += [lo, hi]; bub += [np.zeros(M), np.zeros(M)]; Aeq.append(sm); beq.append([0.0])
d0 = np.zeros((1, nv)); d0[0, off+oid[0]] = 1.0; d0[0, 0] = -3.0; Aeq.append(d0); beq.append([-2.0*N])
Aeq = np.vstack(Aeq); beq = np.concatenate([np.atleast_1d(b) for b in beq]); Aub = np.vstack(Aub); bub = np.concatenate(bub)
c = np.zeros(nv); c[0] = -1.0/N
r = linprog(c, A_ub=csr_matrix(Aub), b_ub=bub, A_eq=csr_matrix(Aeq), b_eq=beq, bounds=(0, None), method='highs')
print(f"N={N} M={M} mode={mode} orbits={no} rows3={len(t3)} status={r.status} maxM2={(-r.fun if r.status==0 else float('nan')):.6f}", flush=True)
