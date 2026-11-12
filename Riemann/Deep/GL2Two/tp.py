"""Two-piece (BCY / KRZ) Levinson-Conrey functional, float64, zeta-normalised variables.
c = c1 + 2 c12 + c2   (BCY 1002.4127 Thm thm:sharp l.264-268; KRZ 1609.03738 Thm unsmoothed l.331-343, l.1745-1748)
  c1  : BCY eq:c1 (l.344-349) = Conrey's c(P1,Q,R,th1)
  c12 : BCY eq:c12 (l.358-366) = KRZ cellell1formula (l.415-421), l = 1
  c2  : BCY eq:c2 (l.374-389)  = KRZ cellellformula (l.392-401), l = 2
GL(2): R = 2 R_KRZ, th_i = nu_i / 2 (KRZ l.1745-1748; Bernard / AT21 eq. (limit-cPQR) for c1).
kappa >= 1 - log(c)/R.
P2 enters only through g = P2'' (P2(0)=P2'(0)=P2''(0)=0); g is expanded in a basis phi_k.
"""
import numpy as np, resource, sys
from numpy.polynomial import chebyshev as C, legendre as Lg

def memguard(cap_mb=280):
    """self-abort above the 300 MB process cap (ru_maxrss is in bytes on macOS)."""
    m = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20
    if m > cap_mb:
        print(f"MEMGUARD: peak RSS {m:.0f} MB > {cap_mb} MB, abort", flush=True); sys.exit(3)

def gl01(n):
    x, w = Lg.leggauss(n)
    return (x + 1) / 2, w / 2

class QCheb:
    """Q(y) = sum_k cz[k] T_k(1-2y)."""
    def __init__(self, cz):
        self.cz = np.asarray(cz, float)
        self.d = [self.cz]
        for k in range(1, 6):
            self.d.append(C.chebder(self.d[-1]))
    def __call__(self, y, k=0):
        return (-2.0) ** k * C.chebval(1 - 2 * np.asarray(y), self.d[k])
    @staticmethod
    def from_h(h):          # Degree2Mol basis: Q = 1 + sum h_n (T_{2n-1}(1-2y) - 1)
        h = np.asarray(h, float); cz = np.zeros(2 * len(h))
        cz[0] = 1 - h.sum()
        for n, hn in enumerate(h, start=1): cz[2 * n - 1] = hn
        return QCheb(cz)
    @staticmethod
    def from_z_mono(a):     # Q = sum_k a_k (1-2y)^k
        return QCheb(C.poly2cheb(np.asarray(a, float)))

class PMono:
    """P(x) = sum_k a[k] x^k (low degree)."""
    def __init__(self, a):
        self.a = np.asarray(a, float)
        self.d = [self.a]
        for k in range(1, 6):
            self.d.append(np.polynomial.polynomial.polyder(self.d[-1]) if len(self.d[-1]) > 1 else np.zeros(1))
    def __call__(self, x, k=0):
        return np.polynomial.polynomial.polyval(np.asarray(x), self.d[k])

def P_sinh(lam, K=12):
    from math import factorial
    a = np.zeros(2 * K + 2)
    for j in range(1, 2 * K + 2, 2): a[j] = lam ** j / factorial(j)
    return PMono(a / a.sum())

def c1(P1, Q, R, th1, n=120):
    v, w = gl01(n)
    q, dq = Q(v), Q(v, 1); U = dq + R * q; e = np.exp(2 * R * v)
    I00 = np.sum(w * e * q * q); I01 = np.sum(w * e * q * U); I11 = np.sum(w * e * U * U)
    u, wu = gl01(60)
    p, dp = P1(u), P1(u, 1)
    a1 = np.sum(wu * dp * dp); a2 = np.sum(wu * p * dp); a3 = np.sum(wu * p * p)
    return 1 + (I00 * a1 + 2 * th1 * I01 * a2 + th1 * th1 * I11 * a3) / th1

# ---- basis for g = P2'' ----
class GBasis:
    """phi_k(s) = s * T_{k-1}(2s-1), k = 1..K  (g(0) = 0); evaluated by chebval (stable)."""
    def __init__(self, K):
        self.K = K
    def ev(self, s, k=None, der=0):
        """array (K,)+s.shape of phi_k^{(der)}(s), der <= 2."""
        s = np.asarray(s, float); z = 2 * s - 1; out = []
        for k in range(1, self.K + 1):
            e = np.zeros(k); e[-1] = 1
            T0 = C.chebval(z, e)
            T1 = 2 * C.chebval(z, C.chebder(e)) if k > 1 else np.zeros_like(z)
            T2 = 4 * C.chebval(z, C.chebder(e, 2)) if k > 2 else np.zeros_like(z)
            out.append([s * T0, T0 + s * T1, 2 * T1 + s * T2][der])
        return np.array(out)
    def P2_coeffs(self, g):
        """monomial coefficients of P2 with P2'' = sum g_k phi_k, P2(0)=P2'(0)=0 (small K only)."""
        tot = np.zeros(self.K + 3)
        for k, gk in enumerate(g, start=1):
            e = np.zeros(k); e[-1] = 1
            pz = np.polynomial.Polynomial(C.cheb2poly(e))
            ps = pz(np.polynomial.Polynomial([-1, 2])) * np.polynomial.Polynomial([0, 1])
            c2 = np.polynomial.polynomial.polyint(ps.coef, 2)
            tot[:len(c2)] += gk * c2
        return tot

def gvec_from_P2(P2coef, K):
    """least-squares (exact if representable) expansion of g = P2'' in GBasis(K)."""
    B = GBasis(K); s, w = gl01(60)
    g = np.polynomial.polynomial.polyval(s, np.polynomial.polynomial.polyder(np.asarray(P2coef, float), 2))
    Phi = B.ev(s, None)
    G = (Phi * w) @ Phi.T; b = (Phi * w) @ g
    return np.linalg.solve(G, b)

# ---- c12 linear functional ----
def c12_vec(P1, Q, R, th1, th2, B, n=40):
    """L_k with c12 = sum_k g_k L_k (BCY eq:c12)."""
    p, wp = gl01(n); q, wq = gl01(n); u, wu = gl01(n)
    PP, QQ, UU = np.meshgrid(p, q, u, indexing='ij')
    W = (wp[:, None, None] * wq[None, :, None] * wu[None, None, :]) * (1 - PP)   # simplex Jacobian
    a = PP; b = (1 - PP) * QQ; s = (1 - PP) * (1 - QQ) * UU
    al = a * UU * th2; be = 1 - b * UU * th2
    A0 = Q(al); A1 = -th1 * (R * Q(al) + Q(al, 1))
    B0 = Q(be); B1 = th1 * (R * Q(be) + Q(be, 1))
    z = 1 - (1 - UU) * th2 / th1
    core = A1 * B1 * P1(z) + (A1 * B0 + A0 * B1) * P1(z, 1) + A0 * B0 * P1(z, 2)
    f = W * UU ** 2 * (1 - UU) * np.exp(R * UU * th2 * (a - b)) * core
    Phi = B.ev(s, None)                     # (K, n, n, n)
    L = np.tensordot(Phi, f, axes=([1, 2, 3], [0, 1, 2]))
    return 4 * (th2 / th1) ** 2 * np.exp(R) * L

# ---- c2 quadratic form ----
def jmul(A, Bj):
    """product of truncated bivariate jets [..., 3, 3]."""
    out = np.zeros(np.broadcast_shapes(A.shape, Bj.shape))
    for i in range(3):
        for j in range(3):
            acc = 0
            for i2 in range(i + 1):
                for j2 in range(j + 1):
                    acc = acc + A[..., i2, j2] * Bj[..., i - i2, j - j2]
            out[..., i, j] = acc
    return out

def qjet(Q, al, be, ga):
    """jet of Q(al + be x + ga y)."""
    fac = [1, 1, 2]
    J = np.zeros(al.shape + (3, 3))
    ders = [Q(al, k) for k in range(5)]
    for i in range(3):
        for j in range(3):
            J[..., i, j] = ders[i + j] * be ** i * ga ** j / (fac[i] * fac[j])
    return J

def c2_mat(Q, R, th2, B, nt=80, nr=24, nu=24, tchunk=40):
    t, wt = gl01(nt); r, wr = gl01(nr); u, wu = gl01(nu)
    K = B.K
    What = np.zeros((nr, nu, nu, 3, 3))
    for ir in range(nr):
        rr = r[ir]
        for t0 in range(0, nt, tchunk):
            T, U, V = np.meshgrid(t[t0:t0 + tchunk], u, u, indexing='ij')         # (nt_c, nu, nv)
            S0 = -rr * (U + V)
            kap = th2 * R * (2 * T - 1)
            E = np.zeros(T.shape + (3, 3))
            e0 = np.exp(2 * R * T + kap * S0)
            cx = kap * (1 - U); cy = kap * (1 - V)
            ex = [np.ones_like(cx), cx, cx * cx / 2]; ey = [np.ones_like(cy), cy, cy * cy / 2]
            for i in range(3):
                for j in range(3):
                    E[..., i, j] = e0 * ex[i] * ey[j]
            Lin = np.zeros(T.shape + (3, 3))
            Lin[..., 0, 0] = 1 / th2 + S0; Lin[..., 1, 0] = 1 - U; Lin[..., 0, 1] = 1 - V
            E = jmul(Lin, E)
            al1 = th2 * U * rr + T * (1 + th2 * S0); be1 = th2 * U + T * th2 * (1 - U); ga1 = -th2 + T * th2 * (1 - V)
            al2 = th2 * V * rr + T * (1 + th2 * S0); be2 = -th2 + T * th2 * (1 - U); ga2 = th2 * V + T * th2 * (1 - V)
            Wj = jmul(jmul(E, qjet(Q, al1, be1, ga1)), qjet(Q, al2, be2, ga2))
            What[ir] += np.tensordot(wt[t0:t0 + tchunk], Wj, axes=(0, 0))
            memguard()
    # X_k[m](r,u): jets of (x+r) phi_k((1-u)(x+r))
    RR, UU = np.meshgrid(r, u, indexing='ij'); c = 1 - UU; s0 = c * RR
    f0 = B.ev(s0, None, 0); f1 = B.ev(s0, None, 1); f2 = B.ev(s0, None, 2)
    X = np.stack([RR * f0, f0 + c * RR * f1, (2 * c * f1 + c * c * RR * f2) / 2], axis=1)   # (K, 3, nr, nu)
    wgt = wr * (1 - r) ** 4
    M = np.zeros((K, K))
    for i in range(3):
        for j in range(3):
            # sum_{r,u,v} wgt[r] wu[u] wv[v] What[r,u,v,i,j] X[k,2-i,r,u] X[l,2-j,r,v]
            Wij = What[..., i, j] * wgt[:, None, None] * wu[None, :, None] * wu[None, None, :]
            M += np.einsum('kru,ruv,lrv->kl', X[:, 2 - i], Wij, X[:, 2 - j], optimize=True)
    M = 8.0 / 3.0 * M
    return (M + M.T) / 2, np.abs(M - M.T).max()
