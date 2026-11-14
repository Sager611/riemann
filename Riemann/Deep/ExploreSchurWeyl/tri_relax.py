#!/usr/bin/env python3
"""Scalable 3-point relaxation (explore50): pair measure nu2 on Z_M and all-index triple measure nu3 on Z_M^2 of an N-multiset
with multiplicities in {1,2}.  max nu2(0)/N = E[sum m^2]/N subject to
 nu2, nu3 >= 0; nu2hat(0)=N^2, nu2hat(k)=k (1<=k<=N), nu2hat(k)>=0 (k>N);
 marginal sum_y nu3(x,y) = N nu2(x);
 coincidence lines L in {x=0, y=0, x=y}: nu2(t) <= nu3|_L(t) <= 2 nu2(t), sum_t nu3|_L(t) = N nu2(0), nu3(0,0) = 3 nu2(0) - 2N;
 Fourier of nu3 (frequencies (a+b,-a,-b), DS range (|a|+|b|+|a+b|)/2 <= N): Haar values (0 if a,b,a+b != 0; N|k| if one is 0; N^3).
 mode 'pair' drops all nu3 Fourier rows except the trivial ones (so it is the pair relaxation + consistent nu3).
usage: tri_relax.py N M [mode]"""
import sys, numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack, hstack, identity
N = int(sys.argv[1]); M = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'tri'
X, Y = np.meshgrid(np.arange(M), np.arange(M), indexing='ij'); X = X.ravel(); Y = Y.ravel(); n3 = M*M
ct = np.cos(2*np.pi*np.arange(M)/M)
fr = np.arange(-(M//2) + 1, M//2 + 1)
rows3, t3 = [], []
for a in fr:
    for b in fr:
        if (abs(a) + abs(b) + abs(a+b)) / 2 > N: continue
        if mode == 'pair' and not (a == 0 and b == 0): continue
        if a == 0 and b == 0: t = N**3
        elif a == 0: t = N*abs(b)
        elif b == 0: t = N*abs(a)
        elif a + b == 0: t = N*abs(a)
        else: t = 0.0
        rows3.append(ct[(a*X + b*Y) % M]); t3.append(t)
A3 = np.array(rows3); t3 = np.array(t3)
# nu2 Fourier
ks = np.arange(0, M//2 + 1); C2 = ct[(np.outer(ks, np.arange(M))) % M]
# variable order: [nu2 (M), nu3 (M^2)]
Z2 = np.zeros((len(t3), M))
Aeq = [np.hstack([Z2, A3])]; beq = [t3]
Aeq.append(np.hstack([C2[:N+1], np.zeros((N+1, n3))])); beq.append(np.concatenate([[N**2], np.arange(1, N+1)]))
Aub = [np.hstack([-C2[N+1:], np.zeros((len(ks)-N-1, n3))])]; bub = [np.zeros(len(ks)-N-1)]
# marginal sum_y nu3(x,y) = N nu2(x)
Mg = np.zeros((M, M + n3))
for x in range(M):
    Mg[x, x] = -N; Mg[x, M + x*M + np.arange(M)] = 1.0
Aeq.append(Mg); beq.append(np.zeros(M))
# coincidence lines
lines = [lambda t: (0, t), lambda t: (t, 0), lambda t: (t, t)]
for L in lines:
    idx = np.array([M + (L(t)[0] % M)*M + (L(t)[1] % M) for t in range(M)])
    lo = np.zeros((M, M + n3)); hi = np.zeros((M, M + n3)); sm = np.zeros((1, M + n3))
    for t in range(M):
        lo[t, t] = 1.0; lo[t, idx[t]] = -1.0          # nu2(t) - line(t) <= 0
        hi[t, idx[t]] = 1.0; hi[t, t] = -2.0          # line(t) - 2 nu2(t) <= 0
    sm[0, idx] = 1.0; sm[0, 0] = -N                   # sum line = N nu2(0)
    Aub += [lo, hi]; bub += [np.zeros(M), np.zeros(M)]; Aeq.append(sm); beq.append([0.0])
d0 = np.zeros((1, M + n3)); d0[0, M] = 1.0; d0[0, 0] = -3.0      # nu3(0,0) - 3 nu2(0) = -2N
Aeq.append(d0); beq.append([-2.0*N])
Aeq = np.vstack(Aeq); beq = np.concatenate([np.atleast_1d(b) for b in beq]); Aub = np.vstack(Aub); bub = np.concatenate(bub)
c = np.zeros(M + n3); c[0] = -1.0/N
r = linprog(c, A_ub=csr_matrix(Aub), b_ub=bub, A_eq=csr_matrix(Aeq), b_eq=beq, bounds=(0, None), method='highs')
print(f"N={N} M={M} mode={mode} rows3={len(t3)} vars={M+n3} status={r.status} maxM2={(-r.fun if r.status==0 else float('nan')):.6f}", flush=True)
