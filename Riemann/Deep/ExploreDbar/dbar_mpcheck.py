# explore59 dbar_mpcheck.py : 40-digit cross-check of U for a split pair (eta=1e-4), K=30, L=20 vs double precision
import mpmath as mp, numpy as np
from dbar_U import U, zz
mp.mp.dps = 40
a = mp.mpf(1) / 2
def K(z, w, L):
    wb = mp.conj(w); D = z - wb
    ph = mp.exp(-1j * L * D)
    I = L if abs(D) < mp.mpf(10) ** -30 else (1 - ph) / (1j * D)
    return ((a * a + wb * wb) * I + (a - 1j * wb) + (a + 1j * wb) * ph) / (2 * a)
g = [mp.mpf(x) for x in zz[:30]]
L = mp.mpf(20); eta = mp.mpf('1e-4')
gg = [g[0] - eta / 2, g[0] + eta / 2] + g[1:]
Z = gg + [-x for x in gg]
n = len(Z)
G = mp.matrix(n, n)
for i in range(n):
    for j in range(n):
        G[i, j] = K(Z[i], Z[j], L)
one = mp.matrix([1] * n)
c = mp.lu_solve(G, one)
u = mp.re(sum(c))
Zd = np.array([float(x) for x in Z])
ud = U(Zd, 20.0)
print('U(mp,40 digits) =', mp.nstr(u, 20), '  U(double) =', repr(ud), '  rel diff =', mp.nstr(abs(u - ud) / u, 3))
