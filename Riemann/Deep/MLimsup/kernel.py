#!/usr/bin/env python3
"""Optimal kernel at a fixed point y (mlimsup).

Maximise  sum_{i<=14400} a_i k(gamma_i) - sum_{14400<i<=N(T1)} 2|c_i| |k(gamma_i)|,
a_i = 2|c_i| cos(gamma_i y + arg c_i) (sign flipped for the liminf point), over k = (g*g~)/||g||^2 with g a step
function on [-T1/2, T1/2] with M cells of width h = 13.55 (T1 = M h).  The zeros above gamma_14400 have unknown
phases at Weisz's points (their ordinates are known only to ~25 digits), so they enter with the worst case.
k is the linear interpolation of the autocorrelation R(d)/R(0) of the cell heights b, so with the worst-case terms
written as -2|c_i| k(gamma_i) (valid while k >= 0 there; checked, and certify.py uses |k| anyway) the objective is
b^T Q b / b^T b with Q symmetric Toeplitz, Q_jl = q(|j-l|)/2 (j != l), Q_jj = q(0), q(d) = sum_i a_i max(0, 1-|gamma_i/h-d|).
Optimum = largest eigenvalue of Q.
Usage: python3 kernel.py phases.npz [ext_zeros.npz]   (phases from certify.py --phases; ext zeros from ext_zeros.py)
Writes kernel_plus.txt and kernel_minus.txt next to this script."""
import sys, os
import numpy as np
from scipy.linalg import toeplitz, eigh
HERE = os.path.dirname(os.path.abspath(__file__))
P = np.load(sys.argv[1])
gam = P['gam']; c = P['c_re'] + 1j * P['c_im']
if len(sys.argv) > 2:
    E = np.load(sys.argv[2]); ge, ce = E['gam'], E['cup']; T1max = float(E['T1max'])
else:
    ge, ce = np.zeros(0), np.zeros(0); T1max = 13550.0
h_num, h_den = 1355, 100
h = h_num / h_den
u = gam / gam[-1]; f = (1 - u) * np.cos(np.pi * u) + np.sin(np.pi * u) / np.pi
Ms = sorted(set([1000] + [int(v) for v in np.arange(1050, int(T1max / h) + 1, 25)] + [int(T1max / h)]))
for w, sg in (('plus', 1.0), ('minus', -1.0)):
    a = sg * 2 * np.abs(c) * np.cos(P['theta_' + w] + np.angle(c))
    best = None
    for M in Ms:
        T1 = M * h
        m_e = ge < T1 + 1          # ext zeros below T1 (+ margin; k = 0 beyond T1 anyway)
        G = np.concatenate([gam, ge[m_e]]); A = np.concatenate([a, -2 * ce[m_e]])
        x = G / h; i0 = np.floor(x).astype(int); fr = x - i0; m = i0 < M
        q = np.zeros(M + 2)
        np.add.at(q, i0[m], (A * (1 - fr))[m]); np.add.at(q, i0[m] + 1, (A * fr)[m])
        col = q[:M] / 2; col[0] = q[0]
        Q = toeplitz(col)
        ev, V = eigh(Q, subset_by_index=[M - 1, M - 1])
        b = V[:, 0] * np.sign(V[:, 0].sum())
        bi = np.rint(b / np.abs(b).max() * 1e6).astype(np.int64)
        bf = bi.astype(float)
        r = np.correlate(bf, bf, 'full')[M - 1:]; r = np.append(r, [0.0, 0.0]) / r[0]
        k = np.where(m, r[np.minimum(i0, M)] * (1 - fr) + r[np.minimum(i0 + 1, M + 1)] * fr, 0.0)
        val = (a * k[:len(gam)]).sum() - (2 * ce[m_e] * np.abs(k[len(gam):])).sum()
        print('%s M=%d T1=%.2f ext zeros=%d  lambda=%.6f  bound=%.6f  min k(ext)=%.2e' % (w, M, T1, m_e.sum(), sg * ev[0], sg * val, k[len(gam):].min() if m_e.any() else 0), flush=True)
        if best is None or val > best[0]:
            best = (val, M, bi)
    val, M, bi = best
    print('%s: JP %.9f  chosen M=%d (T1=%.2f) bound %.9f' % (w, sg * (a * f).sum(), M, M * h, sg * val))
    with open(os.path.join(HERE, 'kernel_%s.txt' % w), 'w') as fh:
        fh.write('# mlimsup kernel for the %s point: g = sum_j b_j * indicator[-T1/2 + j h, -T1/2 + (j+1) h), k = g*g~/||g||^2\n' % w)
        fh.write('# M = %d\n# h_num = %d\n# h_den = %d\n# T1 = M*h = %s\n' % (M, h_num, h_den, M * h_num / h_den))
        fh.write('\n'.join(str(v) for v in bi) + '\n')
