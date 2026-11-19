# fintest.py: falsification test of the finite inequalities behind (H) (isoaudit). numpy, double precision.
# (F0) |I_1| >= 2N - sum K(z-s)^2 - lam|C| - sum_{x in I_2S} min_witness H_tau(y)
# (F)  (1-lam)|I_1| >= (2-lam)N - sum K(z-s)^2 - sum_{distinct pairs} (H_tau(y) - 3 lam)^+
# for every finite conjugation-invariant multiset; lam = lam_C(tau). Run: uv run --no-project --with numpy python3 fintest.py
import numpy as np, sys, time
rng = np.random.default_rng(20261119)
def kern(w):
    A = w / (4*np.sin(w/2))
    def K(xi):   # K(xi) = int f_w(x) e^{-2 pi i xi x} dx, entire; closed form via sinc, complex xi allowed
        xi = np.asarray(xi, dtype=complex)
        a, b = (w + 2*np.pi*xi)/2, (w - 2*np.pi*xi)/2
        sa = np.where(np.abs(a) < 1e-12, 1.0, np.sin(a)/np.where(a == 0, 1, a))
        sb = np.where(np.abs(b) < 1e-12, 1.0, np.sin(b)/np.where(b == 0, 1, b))
        return A*(sa + sb)
    return K
def quadK(w, xi, n=4000):
    x = (np.arange(n) + 0.5)/n - 0.5
    f = np.cos(w*x)/(2*np.sin(w/2)/w)
    return np.sum(f*np.exp(-2j*np.pi*xi*x))/n
def setup(tau, w):
    K = kern(w)
    k0, k1 = K(2*tau).real, K(4*tau).real
    lam = max(2*(1-k0), (6 - 16*k0/3 - 8*k1/3)/3, 0.0)
    thb = np.pi**2*tau**2/2
    assert 4*tau <= 1 and 2*K(tau).real**2 >= 1 and thb < 1 and lam < 1
    H = lambda y: 1 - (1-thb)**2*4*K(1j*y).real**2/(1 + K(2j*y).real)
    return K, lam, H
def check(pts, tau, K, lam, H):
    """pts: list of (complex z, multiplicity); for non-real z only z with Im>0 listed (conjugate added)."""
    Z, M = [], []
    for z, m in pts:
        Z.append(z); M.append(m)
        if abs(z.imag) > 0: Z.append(np.conj(z)); M.append(m)
    Z = np.array(Z); M = np.array(M, dtype=float)
    N = M.sum()
    E = np.sum(np.outer(M, M)*(K(Z[:, None] - Z[None, :])**2)).real
    real = [i for i in range(len(Z)) if Z[i].imag == 0]
    R1 = [i for i in real if M[i] == 1]
    R2 = [i for i in real if M[i] > 1]
    S = [i for i in range(len(Z)) if Z[i].imag > 0]
    I = [i for i in R1 if all(abs(Z[i].real - Z[j].real) > 2*tau for j in R1 if j != i)]
    C = len(R1) - len(I)
    I1, Hs = 0, 0.0
    for i in I:
        wr = [j for j in R2 if abs(Z[i].real - Z[j].real) <= tau]
        ws = [j for j in S if abs(Z[i].real - Z[j].real) <= tau]
        if not wr and not ws: I1 += 1
        elif not wr: Hs += min(H(Z[j].imag) for j in ws)
    F0 = I1 - (2*N - E - lam*C - Hs)
    hsum = sum(max(H(Z[j].imag) - 3*lam, 0.0) for j in S)
    F = (1-lam)*I1 - ((2-lam)*N - E - hsum)
    return F0, F, N
