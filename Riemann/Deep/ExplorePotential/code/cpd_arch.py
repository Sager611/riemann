# cpd_arch.py (explore35, N1 measurement): is Weil's archimedean form conditionally positive (mass-zero f) on wider windows?
# Even sector, CCM cosine basis phi_j (weil_matrix.py, copied from WeilSimple): j >= 1 spans the mass-zero even functions.
# Prints lambda_min of Arch on all even f, on mass-zero even f, of Arch + pole, and of Q - pole on the cosh-orthogonal complement.
import sys, numpy as np
from flint import arb, ctx, fmpq
from weil_matrix import WeilMat
from scipy.optimize import brentq
ctx.prec = 128
def parts(cq, N, sector='even'):
    s = WeilMat(cq); idx = list(range(0, N+1)) if sector == 'even' else list(range(1, N+1))
    n = len(idx); A = np.zeros((n, n)); R = np.zeros((n, n))
    pv = np.array([float(s.pole_vec(j, sector).mid()) for j in idx]); Sv = {j: s.S(j) for j in idx}
    for a_, i in enumerate(idx):
        for b_ in range(a_, n):
            j = idx[b_]
            if i == j: ar = s.arch_diag(i, sector)
            else:
                wi, wj = s.w(i), s.w(j); sg = 1 if (i + j) % 2 == 0 else -1
                ar = -2*s.Nj(i)*s.Nj(j)*sg*(wj*Sv[j] - wi*Sv[i])/(wi*wi - wj*wj) if sector == 'even' else -2*(2/s.L)*sg*(wi*Sv[j] - wj*Sv[i])/(wi*wi - wj*wj)
            pr = arb(0)
            for _, ln, wt in s.atoms: pr += wt*s.C(i, j, ln, sector)
            A[a_, b_] = A[b_, a_] = float(ar.mid()); R[a_, b_] = R[b_, a_] = float((-2*pr).mid())
    return A, R, pv
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
lm = lambda M: float(np.linalg.eigvalsh(M)[0])
def vals(c):
    cq = fmpq(int(round(c*10**6)), 10**6); A, R, pv = parts(cq, N)
    P = 2*np.outer(pv, pv)
    Qp = A + R                                   # Q - pole
    u = pv/np.linalg.norm(pv); Pr = np.eye(len(u)) - np.outer(u, u)
    w, V = np.linalg.eigh(Pr @ Qp @ Pr); w = np.sort(w)
    return lm(A), lm(A[1:, 1:]), lm(A + P), lm(A + P + R), w[1]   # w[0] = 0 (the removed direction)
for c in ((1.1, 1.239, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0) if len(sys.argv) < 3 else (2.5, 3.0)):
    v = vals(c)
    print('c=%.3f a=%.4f | Arch even %.4f | Arch mass-zero %.4f | Arch+pole %.5f | Q %.3e | (Q-pole) on cosh-perp %.3e' % ((c, np.log(c)/2) + v), flush=True)
cz = brentq(lambda c: vals(c)[1], 1.3, 30.0, xtol=1e-5)
print('mass-zero Arch threshold: c = %.5f, a = %.5f (N = %d)' % (cz, np.log(cz)/2, N))
