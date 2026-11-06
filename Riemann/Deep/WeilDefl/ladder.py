# ladder.py (weildefl): is explore17's Phi-ladder the deflated subspace of the certificate matrix T?  NOT part of the proof.
# Reads the dump of cert.py (WD_DUMP): Yt (m x k, orthonormal-ish basis of T's near-null invariant subspace) and the negative columns of B
# (phi_j scaled by sqrt(deff - dl_j); e_k = cos(pi k x/L) scaled by sqrt(-w_k); even sector only here). For v_i = Phi^{(2i)}|_W the vector
# y_i = B^T v_i is T's near-null vector if v_i were M's near-null function (c = B^T psi, DERIVATION.md sec. 1). Prints the cosines of the
# principal angles between span{y_0..y_{K-1}} and span(Yt) for K = k and K = k + 4.
# <v, phi_j>, <v, e_k>: Gauss-Legendre (mpmath) for low frequency, the boundary expansion (integration by parts) for r > 25 * 2 pi c.
import sys, numpy as np, mpmath as mp
mp.mp.dps = 50
d = np.load(sys.argv[1]); c = float(sys.argv[2]); Kx = int(sys.argv[3]) if len(sys.argv) > 3 else 4
Yt = d['Yt']; ctype = d['ctype']; cidx = d['cidx']; cs = d['cscale']; L = mp.log(c); a = L/2
m, k = Yt.shape
NT = 8
def phid(u, J):     # Phi^{(j)}(u), j = 0..J, u >= 0 (z = e^{2u}; recursion (2 z d/dz) on z^a e^{-bz})
    z = mp.e**(2*u); out = [mp.mpf(0)]*(J + 1)
    for n in range(1, NT + 1):
        b = mp.pi*n*n
        for (cf, aa) in ((2*mp.pi**2*n**4, mp.mpf(9)/4), (-3*mp.pi*n*n, mp.mpf(5)/4)):
            poly = [cf]
            for j in range(J + 1):
                out[j] += mp.e**(-b*z)*z**aa*mp.fsum(poly[t]*z**t for t in range(len(poly)))
                new = [mp.mpf(0)]*(len(poly) + 1)
                for t, ct in enumerate(poly):
                    new[t] += 2*(aa + t)*ct; new[t + 1] += -2*b*ct
                poly = new
    return out
K = k + Kx; J = 2*(K - 1) + 14
xg, wg = np.polynomial.legendre.leggauss(400)
X = [(mp.mpf(float(x)) + 1)*a/2 for x in xg]; Wq = [mp.mpf(float(w))*a/2 for w in wg]
D = [phid(x, 2*(K - 1)) for x in X]           # interior values on [0, a]
jet = phid(a, J)                              # boundary jets at L/2
rmax = 25*2*mp.pi*c
def ip(kind, idx, i):    # <v_i, b>, v_i = Phi^{(2i)} on W (even), b = phi_j (normalised) or e_k = cos(pi k x/L)
    if kind == 'phi':
        w = 2*mp.pi*idx/L; Nj = mp.sqrt(1/L) if idx == 0 else mp.sqrt(2/L)
        if w < rmax: return Nj*2*mp.fsum(Wq[q]*D[q][2*i]*mp.cos(w*X[q]) for q in range(len(X)))
        sg = (-1)**idx; return Nj*2*sg*mp.fsum((-1)**t*jet[2*i + 2*t + 1]/w**(2*t + 2) for t in range(6))
    r = mp.pi*idx/L
    if r < rmax: return 2*mp.fsum(Wq[q]*D[q][2*i]*mp.cos(r*X[q]) for q in range(len(X)))
    sk = 1 if ((idx - 1)//2) % 2 == 0 else -1
    return 2*sk*mp.fsum((-1)**t*jet[2*i + 2*t]/r**(2*t + 1) for t in range(6))
Ymat = np.array([[float(cs[col]*ip(ctype[col], int(cidx[col]), i)) for i in range(K)] for col in range(m)])
Qa, _ = np.linalg.qr(Yt)
for KK in (k - 4, k, K):
    Ql, sv = np.linalg.qr(Ymat[:, :KK])[0], None
    s_ = np.linalg.svd(Qa.T @ Ql, compute_uv=False)
    print('c=%g m=%d k=%d  ladder rungs K=%d: cosines of principal angles with span(Yt): %s' % (c, m, k, KK, ' '.join('%.6f' % x for x in s_)), flush=True)
