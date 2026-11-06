# explore20 T1: zeros of the one-sided transform Xi+(z) = 2 int_0^inf K(u) e^{izu} du (Gardner / continuous Enestrom-Kakeya side)
# K = Phi (zeta) or PhiDH; Xi(z) = Xi+(z) + Xi+(-z).  Theorem (EK): K > 0 decreasing on (0, inf) => Xi+ != 0 on Im z >= 0.
import numpy as np, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from offline_np import Phi_np
k5 = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1)
aDH = {0: 0.0, 1: 1.0, 2: k5, 3: -k5, 4: -1.0}
def PhiDH_np(w, N=30):
    w = np.asarray(w, dtype=float); x = np.exp(2*w); s = np.zeros_like(w)
    for n in range(1, N+1):
        if n % 5: s += aDH[n % 5]*n*np.exp(1.5*w - np.pi*n*n*x/5)
    return s
xg, wg = np.polynomial.legendre.leggauss(24)
edges = np.linspace(0, 4.5, 121)
NODES = (((edges[:-1] + edges[1:])/2)[:, None] + (np.diff(edges)/2)[:, None]*xg[None, :]).ravel()
WTS = ((np.diff(edges)/2)[:, None]*wg[None, :]).ravel()
def make_onesided(K):
    KV = K(NODES)*WTS
    def F(z, k=0, chunk=256):
        z = np.atleast_1d(np.asarray(z, dtype=complex)); out = np.empty(z.shape, dtype=complex)
        w = KV*(1j*NODES)**k
        for i in range(0, z.size, chunk):
            zz = z[i:i+chunk]
            out[i:i+chunk] = 2*(np.exp(1j*np.multiply.outer(zz, NODES))*w).sum(-1)
        return out
    return F
def find_zeros(F, xs, ys):
    X, Y = np.meshgrid(xs, ys); Z = X + 1j*Y
    A = np.abs(F(Z.ravel())).reshape(Z.shape)
    cands = [Z[i, j] for i in range(1, A.shape[0]-1) for j in range(1, A.shape[1]-1) if A[i, j] <= A[i-1:i+2, j-1:j+2].min()]
    zs = []
    for z in cands:
        for it in range(60):
            f = F(z)[0]; df = F(z, 1)[0]
            dz = f/df; z = z - dz
            if abs(dz) < 1e-13*(1 + abs(z)): break
        if abs(F(z)[0]) < 1e-10 and xs[0] - 1 <= z.real <= xs[-1] + 1 and ys[0] - 1 <= z.imag <= ys[-1] + 1:
            if all(abs(z - w) > 1e-6 for w in zs): zs.append(z)
    return sorted(zs, key=lambda w: w.real)
