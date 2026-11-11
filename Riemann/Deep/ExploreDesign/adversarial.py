"""Adversarial test (iii): Weil form of zeta vs a xi-like function with the first two zeros merged into an off-line pair,
gamma_1, gamma_2 -> (gamma_1+gamma_2)/2 +- i delta (zero count kept; xi_m = xi (s-r+)(s-r-)/((s-r1)(s-r2)) and conjugates).
Even/odd sectors of the P = 2L lattice form (Chekmachev's Thm 3.2 basis), Galerkin in phi_j, j < J, odd nodes k <= Kmax (float64).
Q_m(f) = Q(f) - |F(g1)|^2 - |F(g2)|^2 + 2 Re F(gm - i d) conj F(gm + i d)  (+ the mirror zeros at -gamma).
Lattice weight of one zero gamma (complex allowed): int_{-L}^{L} e^{i gamma x} e^{-i r_k x} dx = 2 sin((gamma - r_k) L)/(gamma - r_k)."""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from wsym import mu

def zero_weight(gam, r, L):
    z = gam - r
    out = np.where(np.abs(z) < 1e-12, 2 * L + 0j, 2 * np.sin(z * L) / np.where(np.abs(z) < 1e-12, 1, z))
    return out

def sector_matrix(c, sector, J, Kmax, dmu=None):
    L = math.log(c)
    j = np.arange(0, J) if sector == 'even' else np.arange(1, J + 1)
    k = np.arange(1, Kmax + 1, 2)
    rk = math.pi * k / L
    sk = (-1.0) ** ((k - 1) // 2)
    if sector == 'even':
        Nj = np.where(j == 0, L ** -0.5, (2 / L) ** 0.5)
        Aj = Nj * (L / math.pi) * (-1.0) ** j
        bk = 2 * k * sk
        yv = -L * L / (4 * math.pi ** 2); bv = -(L / math.pi) * math.sinh(L / 4)
    else:
        Nj = np.full(len(j), (2 / L) ** 0.5)
        Aj = Nj * (-1.0) ** j * 4 * j * (L / math.pi)
        bk = sk.astype(float)
        yv = -L * L / (4 * math.pi ** 2); bv = math.sinh(L / 4)
    G = Aj[:, None] * bk[None, :] / (k[None, :] ** 2 - 4.0 * j[:, None] ** 2)       # <phi_j, e_k>
    pv = Aj * bv / (yv - 4.0 * j ** 2)                                              # <phi_j, v>
    me = mu(2 * math.pi * j / L, c); mo = mu(rk, c)
    if dmu is not None:
        me = me + dmu(2 * math.pi * j / L); mo = mo + dmu(rk)
    Qm = np.diag(0.5 * me) + (G * (mo / L)[None, :]) @ G.T
    Qm += (2.0 if sector == 'even' else -2.0) * np.outer(pv, pv)
    return Qm

if __name__ == '__main__':
    import mpmath
    g1 = float(mpmath.zetazero(1).imag); g2 = float(mpmath.zetazero(2).imag)
    gm = 0.5 * (g1 + g2)
    for c in [2, 3, 5, 10]:
        L = math.log(c)
        for delta in [0.3, 0.05]:
            def dmu(r, d=delta):
                tot = np.zeros_like(r, dtype=complex)
                for s in (1, -1):
                    tot -= zero_weight(s * g1, r, L) + zero_weight(s * g2, r, L)
                    tot += zero_weight(s * gm + 1j * d, r, L) + zero_weight(s * gm - 1j * d, r, L)
                return tot.real
            for sector in ['even', 'odd']:
                J, Kmax = 40, 4001
                lz = np.linalg.eigvalsh(sector_matrix(c, sector, J, Kmax))[0]
                lm = np.linalg.eigvalsh(sector_matrix(c, sector, J, Kmax, dmu))[0]
                print(f'c={c} delta={delta} {sector}: lambda_min zeta {lz:.3e}   merged off-line pair {lm:.3e}')
    # the same perturbation seen by the patch region (t >= 3000, c = 47): size of dmu
    L = math.log(47)
    t = np.linspace(3000, 3100, 2001)
    def dmu47(r, d=0.3):
        tot = np.zeros_like(r, dtype=complex)
        for s in (1, -1):
            tot -= zero_weight(s * g1, r, L) + zero_weight(s * g2, r, L)
            tot += zero_weight(s * gm + 1j * d, r, L) + zero_weight(s * gm - 1j * d, r, L)
        return tot.real
    print('c=47: max |Delta mu| on t in [3000, 3100] =', float(np.max(np.abs(dmu47(t)))))
