# Shared machinery for the re-optimised HB 2610.12803 table (p = 61/20). Float model + exact Bernstein algebra.
import math
from fractions import Fraction as F
from math import comb
import numpy as np

P_ = F(61, 20); B_ = P_ * P_
N1 = 65536

# ---------------- float scalars (model only) ----------------
_n = np.arange(2, N1 + 1, dtype=np.float64); _ln = np.log(_n)
_sieve = np.ones(N1 + 1, dtype=bool); _sieve[:2] = False
for _i in range(2, int(N1 ** 0.5) + 1):
    if _sieve[_i]: _sieve[_i * _i::_i] = False
_pp_n, _pp_j, _pp_ll = [], [], []
for _l in np.nonzero(_sieve)[0]:
    _m, _j = int(_l), 1
    while _m <= N1:
        _pp_n.append(_m); _pp_j.append(_j); _pp_ll.append(math.log(int(_l))); _m *= int(_l); _j += 1
_pp_ln = np.log(np.array(_pp_n, dtype=np.float64)); _pp_j = np.array(_pp_j, dtype=np.float64); _pp_ll = np.array(_pp_ll)
_logN1 = math.log(N1)
def Tf(j, v):
    return N1 ** (1 - v) * sum(math.factorial(j) // math.factorial(j - k) * _logN1 ** (j - k) / (v - 1) ** (k + 1) for k in range(j + 1))
pf = 61 / 20; Bf = pf * pf; sigp = 81 / 40; Nst = 69098.0; Lst = 2 * math.log(Nst) + 3
c1f = 0.5 + 1e-10; c2f = 0.25 + 1e-10
def scalars_f(tm, tp):
    u = sigp + 5.57 * tm; u8 = sigp + 2.78 * tm; e = (pf - 1) / 2 - 2.78 * tm
    npw = np.exp(-u * _ln)
    Z = 1 + npw.sum() + Tf(0, u)
    ppw = np.exp(-u * _pp_ln)
    psi = (_pp_ll * ppw).sum() + Tf(1, u); psi2 = (_pp_j * _pp_ll ** 2 * ppw).sum() + Tf(2, u)
    hm1 = np.expm1(tp * _ln * _ln / 4) * npw
    d0 = hm1.sum() + Tf(0, u8); d1 = (hm1 * _ln).sum() + Tf(1, u8); d2 = (hm1 * _ln * _ln).sum() + Tf(2, u8)
    r0 = 1.01 * Nst ** (-pf) * (Nst ** (e + 1) / (e + 1) + Nst ** e)
    return dict(Z=Z, psi=psi, psi2=psi2, d0=d0, d1=d1, d2=d2, r0=r0, r1=Lst * r0, r2=(Lst * Lst + 1e-9) * r0)

def CG(s, one=1.0, c1=None, c2=None, eta=None, k55716=None, e6=None, e9=None, e10=None):
    # generic: works for floats (defaults) and Fractions (pass exact constants)
    if c1 is None:
        c1, c2 = c1f, c2f; eta = (1e-6, 4e-5, 16e-4); k55716 = 5.5716; e6, e9, e10 = 1e-6, 1e-9, 2e-10
    Z, ps, ps2 = s['Z'], s['psi'], s['psi2']
    e0 = s['d0'] + s['r0'] + eta[0]; ep0 = Z * e0
    e1 = c1 * s['d1'] + s['r1'] + eta[1]; ep1 = Z * (e1 + c1 * ps * e0)
    e2 = c2 * s['d2'] + s['r2'] + eta[2] + e9
    ep2 = Z * (e2 + 2 * c1 * ps * e1 + c2 * (ps2 + ps * ps) * e0)
    C = k55716 - e6 - c1 * ps - ep1 / (one - ep0)
    G = c2 * ps2 + ep2 / (one - ep0) + (ep1 / (one - ep0)) ** 2 + e10
    return C, G, ep0

def CG_exact(s):
    return CG(s, one=F(1), c1=F(1, 2) + F(1, 10 ** 10), c2=F(1, 4) + F(1, 10 ** 10),
              eta=(F(1, 10 ** 6), F(4, 10 ** 5), F(16, 10 ** 4)), k55716=F(55716, 10 ** 4),
              e6=F(1, 10 ** 6), e9=F(1, 10 ** 9), e10=F(2, 10 ** 10))

# ---------------- kernel polynomial Phi and compactified Psi ----------------
def _pmul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            out[(i + k, j + l)] = out.get((i + k, j + l), 0) + x * y
    return out
def _padd(*ps):
    out = {}
    for c, a in ps:
        for key, x in a.items(): out[key] = out.get(key, 0) + c * x
    return out

def psi_matrix(A2, lam, mu, B=B_):
    """5x5 coefficients (U^i V^j) of Psi(U,V) = (1-U)^4/B^4 * Phi(BU/(1-U), A2 V). Exact if inputs are Fractions."""
    one = B / B
    Dp = {(1, 0): one}; rp = {(0, 1): one}
    w = _padd((1, Dp), (B, {(0, 0): one}), (-1, rp))
    Dl = _padd((1, _pmul(w, w)), (4, _pmul(Dp, rp)))
    wD = _pmul(w, Dl)
    Ncal = _padd((2 * lam, wD), (-2 * mu, wD), (-4 * mu * B, Dl), (8 * mu * B, _pmul(w, w)))
    fac = _padd((1, Dp), (A2, {(0, 0): one}), (3, rp))
    Phi = _padd((2, _pmul(Dl, Dl)), (-1, _pmul(fac, Ncal)))
    # Psi: D^i r^j -> B^(i-4) A2^j U^i (1-U)^(4-i) V^j
    M = [[0 * one for _ in range(5)] for _ in range(5)]
    for (i, j), c in Phi.items():
        assert i <= 4 and j <= 4
        if c == 0: continue
        cc = c * B ** (i - 4) * A2 ** j if i <= 4 else None
        # (1-U)^(4-i) * U^i
        for k in range(4 - i + 1):
            M[i + k][j] += cc * comb(4 - i, k) * (-1) ** k
    return M

def _subs_mat(m, l, one):
    # S[i][k] = coeff of u^k in ((l+u)/m)^i
    return [[(comb(i, k) * l ** (i - k)) * one / (m ** i) if k <= i else 0 * one for k in range(5)] for i in range(5)]
_W = [[F(comb(a, k), comb(4, k)) if k <= a else F(0) for k in range(5)] for a in range(5)]

def bernstein_coeffs(M, mU, mV, one):
    """All Bernstein coefficients (list) of Psi on the mU x mV grid of sub-rectangles of [0,1]^2."""
    out = []
    W = [[one * x for x in row] for row in _W]
    for l in range(mU):
        S = _subs_mat(mU, l, one)
        # A1[k][j] = sum_i S[i][k] M[i][j]
        A1 = [[sum(S[i][k] * M[i][j] for i in range(5)) for j in range(5)] for k in range(5)]
        for jv in range(mV):
            Sv = _subs_mat(mV, jv, one)
            A = [[sum(A1[k][j] * Sv[j][nn] for j in range(5)) for nn in range(5)] for k in range(5)]
            # beta = W A W^T
            WA = [[sum(W[a][k] * A[k][nn] for k in range(5)) for nn in range(5)] for a in range(5)]
            for a in range(5):
                for b in range(5):
                    out.append(sum(WA[a][nn] * W[b][nn] for nn in range(5)))
    return out

def bern_affine_f(A2, mU, mV):
    b0 = np.array(bernstein_coeffs(psi_matrix(float(A2), 0.0, 0.0, Bf), mU, mV, 1.0))
    bl = np.array(bernstein_coeffs(psi_matrix(float(A2), 1.0, 0.0, Bf), mU, mV, 1.0)) - b0
    bm = np.array(bernstein_coeffs(psi_matrix(float(A2), 0.0, 1.0, Bf), mU, mV, 1.0)) - b0
    return b0, bl, bm

def alpha_gamma(C, G, q, p=pf, B=Bf):
    return C / p - 2 / (B - q), G + 2 * (B + q) / (B - q) ** 2

def lp_best(A2, C, G, mU=16, mV=2, margin=0.0):
    from scipy.optimize import linprog
    b0, bl, bm = bern_affine_f(A2, mU, mV)
    al, ga = alpha_gamma(C, G, A2)
    res = linprog([-al, ga], A_ub=np.stack([-bl, -bm], 1), b_ub=b0 - margin, bounds=[(0, None), (0, None)], method='highs')
    lam, mu = res.x
    return lam, mu, 4 * (al * lam - ga * mu)

def kexact(lam, mu, C, G, q, p=P_, B=B_):
    return 4 * (lam * C / p - mu * G - 2 * lam / (B - q) - 2 * mu * (B + q) / (B - q) ** 2)

def E(x, n):
    s = F(0); term = F(1)
    for j in range(n + 1):
        if j > 0: term = term * x / j
        s += term
    return s
