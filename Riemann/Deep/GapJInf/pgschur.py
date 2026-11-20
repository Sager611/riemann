"""gapjinf: per-mode (Schur) form of the order-J per-gap problem; exact elimination of every sine mode k <= KM.
Coordinates as pgq.py: per mode k (parity par), (beta_k, tau) with tau = boundary data (sig for odd k, del for even k);
per-mode quadratic data A_v(k), B_v(k) (L), C_v(k) (LxL) for each form variable v (a_0..a_J, Phi_ml), plus constants
(the del^2/g terms). Given a form x: A = x.A_v, B = x.B_v, C = x.C_v; min over beta_k gives R = sum_k (C - B B^T/A) + const.
Valid at (g, par)  <=>  A(k) > 0 for all k and R >= 0."""
import numpy as np
from scipy.special import zeta

PI = np.pi


def nj(j):
    return 0 if j == 0 else (j - 1) // 2


def mode_data(g, J, par, KM=20000):
    L = (J - 1) // 2
    ks = np.arange(1, KM + 1); ks = ks[ks % 2 == par]
    w = ks * PI / g
    sgn = 1.0 if par == 1 else -1.0
    c = np.zeros((len(ks), L + 1))
    for l in range(1, L + 1):
        c[:, l] = sgn * (-1) ** l * (2 / g) * w ** (-2.0 * l - 1)
    tkeys = [(m, l) for m in range(0, (J - 2) // 2 + 1) for l in range(1, L + 1)]
    nv = J + 1 + len(tkeys)
    A = np.zeros((nv, len(ks))); B = np.zeros((nv, len(ks), L)); C = np.zeros((nv, len(ks), L, L)); K0 = np.zeros((nv, L, L))
    for j in range(J + 1):
        n = nj(j)
        r = c[:, 1:].copy(); r[:, :n] = 0.0                   # u_{k,n} = beta + sum_{l>n} c_l tau_l
        f = (g / 2) * w ** (2.0 * j)
        A[j] = f; B[j] = f[:, None] * r; C[j] = f[:, None, None] * r[:, :, None] * r[:, None, :]
        if j % 2 == 1 and j >= 3 and par == 0:
            mm = (j - 1) // 2
            if mm <= L:
                K0[j, mm - 1, mm - 1] += 1.0 / g
    for i, (m, l) in enumerate(tkeys):
        v = J + 1 + i
        r = c[:, 1:].copy(); r[:, :m] = 0.0                   # u_{k,m}
        f = -(g / 2) * (-1) ** (m + l) * w ** (2.0 * (m + l + 1))
        e = np.zeros((len(ks), L)); e[:, l - 1] = c[:, l]      # eps_{k,l} as a vector in tau
        B[v] = 0.5 * f[:, None] * e
        Cv = f[:, None, None] * r[:, :, None] * e[:, None, :]
        C[v] = 0.5 * (Cv + np.transpose(Cv, (0, 2, 1)))
        if m >= 1 and par == 0:
            K0[v, m - 1, l - 1] += 0.5 / g; K0[v, l - 1, m - 1] += 0.5 / g
    return dict(ks=ks, w=w, A=A, B=B, C=C, K0=K0, L=L, nv=nv, tkeys=tkeys, g=g, par=par)


def reduce(D, x):
    """returns (min_k A(k), R (LxL), beta-coefficients matrix Bt = -B/A (KxL))"""
    A = x @ D['A']; B = np.tensordot(x, D['B'], axes=(0, 0)); C = np.tensordot(x, D['C'], axes=(0, 0))
    K0 = np.tensordot(x, D['K0'], axes=(0, 0))
    Amin = A.min()
    if D['L'] == 0:
        return Amin, None, None
    Sk = np.abs(x) @ np.abs(D['A'])
    Ap = np.maximum(A, 1e-13 * Sk)
    terms = C - B[:, :, None] * B[:, None, :] / Ap[:, None, None]
    A = Ap
    R = terms.sum(axis=0) + K0
    # tail k > KM: terms ~ const/k^2 -> sum_{k>KM, step 2} ~ last * KM/2
    R = R + terms[-1] * len(D['ks'])
    return Amin, 0.5 * (R + R.T), -B / A[:, None]


def piece_vector(D, tau, Bt):
    """cut vector w_v = sum_k [A_v beta^2 + 2 beta B_v.tau + tau C_v tau] + tau K0_v tau, beta_k = Bt_k . tau"""
    beta = Bt @ tau
    w = D['A'] @ (beta ** 2) + 2 * np.einsum('vkl,k,l->v', D['B'], beta, tau) + np.einsum('vkab,a,b->v', D['C'], tau, tau)
    return w + np.einsum('vab,a,b->v', D['K0'], tau, tau)
