"""gapjinf: amplifier step of the joint loop: for a saved form (logs/form_*.json, a normalised a1 = 1), the best amplifier is the
top eigenvector of (C0, sum_{j>=1} a_j C_2j(v)) (h = -a0 * lambda_max), v optimised by a 1-D search; writes an amp.json in
explore53's convention (A = R^{-1} y * norm, A_0 = 1) so that ampmom.load can read it.
usage: reopt.py forms.npz form.json out_amp.json"""
import sys, json
from math import comb
import numpy as np
from scipy.optimize import minimize_scalar


def shift(C, dv, kap):
    t = -dv * kap
    return [sum(comb(k, i) * C[i] * t ** (k - i) for i in range(k + 1)) for k in range(len(C))]


def lmax(A, B):
    L = np.linalg.cholesky(B); Li = np.linalg.inv(L)
    ev, U = np.linalg.eigh(Li @ A @ Li.T)
    return ev[-1], Li.T @ U[:, -1]


D = np.load(sys.argv[1]); F = json.load(open(sys.argv[2]))
C = list(D['C']); R = D['R']; norm = D['norm']; kap = float(D['kap']); v0 = float(D['v0'])
a = np.array(F['a']); J = F['J']
def h(dv):
    Cs = shift(C[:2 * J + 1], dv, kap)
    return -a[0] * lmax(Cs[0], sum(a[j] * Cs[2 * j] for j in range(1, J + 1)))[0]
rr = minimize_scalar(lambda x: -h(x), bounds=(-0.01, 0.01), method='bounded', options={'xatol': 1e-10})
dv = rr.x; Cs = shift(C[:2 * J + 1], dv, kap)
lam, y = lmax(Cs[0], sum(a[j] * Cs[2 * j] for j in range(1, J + 1)))
co = np.linalg.solve(R, y); A = co * norm; A = A / A[0]
print('kappa=%.4f J=%d: form h on old amplifier %.7f -> best amplifier h = %.7f at v = %.6f' % (kap, J, F['h'], -a[0] * lam, v0 + dv))
json.dump({'v': '%.6f' % (v0 + dv), 'A': ['%.17g' % x for x in A], 'a': ['%.17g' % x for x in a], 'kappa': kap, 'J': J}, open(sys.argv[3], 'w'), indent=1)
