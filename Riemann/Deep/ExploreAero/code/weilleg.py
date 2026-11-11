# weilleg.py: Weil's quadratic form for zeta on even functions on the window [-a, a] (c = e^{2a}) in the orthonormal even
# Legendre basis p_n = sqrt((2n+1)/(2a)) P_n(x/a), n = 0, 2, ..., 2K-2 (double precision; orientation only, not rigorous).
#   Q(f) = 2<f,cosh(x/2)>^2 + (1/2pi) int (Re psi(1/4+ir/2) - log pi)|F(r)|^2 dr - 2 sum_{k<c} Lambda(k) k^{-1/2} C_f(log k)
#   arch in x-space (psi integral rep.): (-gamma - log pi) g(0) + int_0^inf [g(0)e^{-2x} - e^{-x/2} g(x)] 2dx/(1-e^{-2x}),
#   g = C_f = f * f~ (even, supported on [-2a, 2a]); tail int_{2a}^inf g(0) 2e^{-2x}/(1-e^{-2x}) dx = -g(0) log(1-e^{-4a}).
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import spherical_jn, spherical_in
EULER = 0.57721566490153286061; LOGPI = np.log(np.pi)
def vonmangoldt(k):
    for p in range(2, k + 1):
        if k % p == 0:
            m = k
            while m % p == 0: m //= p
            return np.log(p) if m == 1 else 0.0
    return 0.0
class WeilLeg:
    def __init__(s, a, K, nin=None, nout=None):
        s.a = a; s.K = K; s.deg = 2*np.arange(K); s.nrm = np.sqrt((2*s.deg + 1)/(2*a))
        dmax = s.deg[-1]
        s.nin = nin or (dmax + 40); s.nout = nout or (2*dmax + 60)
        s.ti, s.wi = Lg.leggauss(s.nin)
    def vals(s, x):   # K x len(x) basis values (zero outside [-a, a])
        t = np.asarray(x)/s.a
        V = Lg.legvander(np.clip(t, -1, 1), s.deg[-1])[:, s.deg].T*s.nrm[:, None]
        V[:, np.abs(t) > 1] = 0.0
        return V
    def C(s, x):      # C_nm(x) = int p_n(y) p_m(y - x) dy, 0 <= x <= 2a
        lo, hi = x - s.a, s.a
        y = lo + (hi - lo)*(s.ti + 1)/2; w = s.wi*(hi - lo)/2
        A = s.vals(y); B = s.vals(y - x)
        return (A*w) @ B.T
    def pole(s):
        return s.nrm*2*s.a*spherical_in(s.deg, s.a/2)
    def ft(s, z):     # p_n^(z) = int p_n(x) e^{izx} dx = nrm * 2a * i^n j_n(a z)
        z = np.atleast_1d(z)
        J = np.array([spherical_jn(int(n), s.a*z) for n in s.deg])
        return (s.nrm*2*s.a*((-1.0)**(s.deg//2)))[:, None]*J
    def build(s):
        K = s.K; C0 = s.C(0.0)
        to, wo = Lg.leggauss(s.nout); x = s.a*(to + 1); wx = wo*s.a
        A = (-EULER - LOGPI)*C0 - C0*np.log(1 - np.exp(-4*s.a))
        for xi, wi in zip(x, wx):
            A += wi*(C0*np.exp(-2*xi) - np.exp(-xi/2)*s.C(xi))*2/(1 - np.exp(-2*xi))
        P = np.zeros((K, K)); k = 2
        while np.log(k) < 2*s.a:
            L = vonmangoldt(k)
            if L > 0: P += 2*L/np.sqrt(k)*s.C(np.log(k))
            k += 1
        l = s.pole()
        s.A, s.P, s.l, s.G = (A + A.T)/2, (P + P.T)/2, l, C0
        s.Q = 2*np.outer(l, l) + s.A - s.P
        return s
