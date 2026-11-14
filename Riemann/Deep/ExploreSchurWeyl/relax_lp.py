#!/usr/bin/env python3
"""Pair-measure LP relaxation (finite-N analogue of the CGdL/Cohn-Elkies LP), on Z_M.
max nu({0})/N  s.t. nu >= 0 on Z_M (symmetric), nuhat(0)=N^2, nuhat(k)=k (1<=k<=N), nuhat(k) >= 0 (N<k<=M/2).
usage: relax_lp.py N M"""
import sys, numpy as np
from scipy.optimize import linprog
N = int(sys.argv[1]); M = int(sys.argv[2])
xs = np.arange(0, M//2 + 1)                       # symmetric nu: nu_x = nu_{-x}; variables y_x = nu_x + nu_{-x} (x>0), y_0 = nu_0
ks = np.arange(0, M//2 + 1)
Cmat = np.cos(2*np.pi*np.outer(ks, xs)/M)         # nuhat(k) = sum_x y_x cos(2 pi k x/M)
Aeq = Cmat[:N+1]; beq = np.concatenate([[N**2], np.arange(1, N+1)])
Aub = -Cmat[N+1:]; bub = np.zeros(len(ks) - N - 1)
c = np.zeros(len(xs)); c[0] = -1.0/N
r = linprog(c, A_ub=Aub if len(bub) else None, b_ub=bub if len(bub) else None, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
print(f"N={N} M={M} relax maxM2={-r.fun:.6f} status={r.status}")
