# explore42 s14: zeta-free floor F(lam) = dist^2(chi_1, W_lam), W_lam = {f in L^2(0,1): V f = 0 on (0,lam)} (contains C_lam).
# Reduction (DERIVATION.md): F = <q, (I - T)^{-1} q>_{L^2(0,lam)},  q = C J chi_1,  q(v) = -2[(1 + log a) Ci(a) + J(a)], a = 2 pi v,
# J(a) = int_a^inf cos(u) log(u)/u du;  T(v,v') = sin(2pi(v-v'))/(pi(v-v')) + sin(2pi(v+v'))/(pi(v+v')) (even prolate concentration, band (0,1)).
# eps = ||q||^2 = ||(1-Q_lam) V chi_1||^2 (check against s11).
import sys, math
import numpy as np
import mpmath as mp
from scipy.special import sici
mp.mp.dps = 25
J1 = float(mp.quadosc(lambda u: mp.cos(u) * mp.log(u) / u, [1, mp.inf], omega=1))
def Jf(a):
    # int_a^1 cos u log u/u = -log(a)^2/2 + int_a^1 (cos u - 1) log u/u du
    r = float(mp.quad(lambda u: (mp.cos(u) - 1) * mp.log(u) / u, [a, 1]))
    return J1 - math.log(a) ** 2 / 2 + r
def q(v):
    a = 2 * math.pi * v
    ci = sici(a)[1]
    return -2 * ((1 + math.log(a)) * ci + Jf(a))
def Tk(v, w):
    d = v[:, None] - w[None, :]; s = v[:, None] + w[None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        t1 = np.where(np.abs(d) < 1e-14, 2.0, np.sin(2 * np.pi * d) / (np.pi * d))
    t2 = np.sin(2 * np.pi * s) / (np.pi * s)
    return t1 + t2
def F(N, M, p=3):
    # h ranges over H3 = (1-M')^3 L^2(0,lam) = {h ⊥ 1, log(lam/v), log^2(lam/v)} (so that g = (1-M)^3 h stays in L^2(0,lam))
    lam = 1.0 / N
    x, wx = np.polynomial.legendre.leggauss(M); x = (x + 1) / 2; wx = wx / 2
    v = lam * x ** p; w = wx * lam * p * x ** (p - 1)
    qv = np.array([q(t) for t in v])
    K = Tk(v, v)
    sw = np.sqrt(w)
    # work in the w-weighted coordinates z = sqrt(w) h
    Cc = np.stack([sw * np.log(lam / v) ** k for k in range(3)], axis=1)       # constraint vectors in z-coords
    Q, _ = np.linalg.qr(Cc, mode='complete'); Z = Q[:, 3:]                      # orthonormal basis of the constraint complement
    Kz = sw[:, None] * K * sw[None, :]
    Az = np.eye(M) - Kz
    b = Z.T @ (sw * qv)
    Ared = Z.T @ Az @ Z
    Fv = float(b @ np.linalg.solve(Ared, b))
    eps = float(b @ b)                                                           # sup |<q,h>|^2/||h||^2 over H3
    mu = np.linalg.eigvalsh(Z.T @ Kz @ Z)
    return eps, Fv, mu[-3:][::-1]
print('J(1) =', J1)
for N in [2, 3, 4, 5, 6, 8, 10, 20, 40, 50, 100]:
    for M in [40, 80]:
        e, Fv, mu = F(N, M)
        print(f'N={N:4d} M={M:3d} eps={e:.10e} F={Fv:.10e} F-eps={Fv-e:.4e} (F-eps)logN={(Fv-e)*math.log(N):.4e} F/eps={Fv/e:.5f} top mu={np.array2string(mu, precision=4)}', flush=True)
