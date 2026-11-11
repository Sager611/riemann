"""Baluyot's mollifier (thesis (1.2.6)-(1.2.7)) and the diagonal main terms.

W(xi; sigma) = (1/2 pi i) int_(r) Gamma(s+1) zeta(2sigma+s)/zeta(2sigma) zeta(1+s) xi^s ds
            = (1/(xi zeta(2sigma))) sum_b b^{1-2sigma} / (e^{b/xi} - 1)        (termwise inverse Mellin)
            = 1 + A xi^{1-2sigma} + sum_k c_k xi^{-k}                           (residues, xi >= 1)
A = Gamma(2-2sigma) zeta(2-2sigma)/zeta(2sigma),
c_k = (-1)^{k-1}/(k-1)! zeta(1-k) zeta(2sigma-k)/zeta(2sigma).
M(s) = sum_{n <= N} mu(n) W(y/n)/W(y) n^{-s},  N = floor(y^{1+eps}).
M1 = zeta(2sigma) sum_{h,k} a_h a_k / [h,k]^{2sigma}            (thesis (3.3.9) with I = 1)
S2 = sum_{h,k} a_h a_k (h,k)^{2-2sigma}/(hk)                    (inner sum of M2, (3.3.11))
"""
import numpy as np
import mpmath

K_ASYM = 40


def mobius(N):
    mu = np.ones(N + 1, dtype=np.int8)
    mu[0] = 0
    isp = np.ones(N + 1, dtype=bool)
    isp[:2] = False
    for p in range(2, N + 1):
        if isp[p]:
            isp[2 * p::p] = False
            mu[p::p] *= -1
            if p * p <= N:
                mu[p * p::p * p] = 0
    return mu, np.nonzero(isp)[0]


class Wfun:
    def __init__(self, sigma):
        mpmath.mp.dps = 30
        self.sigma = sigma
        z2s = mpmath.zeta(2 * sigma)
        self.z2s = float(z2s)
        self.A = float(mpmath.gamma(2 - 2 * sigma) * mpmath.zeta(2 - 2 * sigma) / z2s)
        self.ck = np.array([float((-1) ** (k - 1) / mpmath.factorial(k - 1) * mpmath.zeta(1 - k)
                                  * mpmath.zeta(2 * sigma - k) / z2s) for k in range(1, K_ASYM + 1)])

    def series(self, xi):
        xi = np.asarray(xi, dtype=float)
        out = np.zeros_like(xi)
        bmax = int(np.ceil(60 * xi.max())) + 10 if xi.size else 1
        for b in range(1, bmax + 1):
            out += b ** (1 - 2 * self.sigma) / np.expm1(b / xi)
        return out / (xi * self.z2s)

    def asym(self, xi):
        xi = np.asarray(xi, dtype=float)
        out = 1 + self.A * xi ** (1 - 2 * self.sigma)
        inv = 1.0 / xi
        p = np.ones_like(xi)
        for c in self.ck:
            p = p * inv
            out = out + c * p
        return out

    def __call__(self, xi):
        xi = np.asarray(xi, dtype=float)
        out = np.empty_like(xi)
        big = xi >= 1.0
        if big.any():
            out[big] = self.asym(xi[big])
        if (~big).any():
            out[~big] = self.series(xi[~big])
        return out


def coeffs(y, sigma, N):
    mu, primes = mobius(N)
    W = Wfun(sigma)
    n = np.arange(1, N + 1, dtype=float)
    a = np.zeros(N + 1)
    nz = mu[1:] != 0
    a[1:][nz] = mu[1:][nz] * W(y / n[nz]) / W(np.array([float(y)]))[0]
    return a, mu, primes, W


def jordan_sqfree(N, alpha, primes):
    J = np.ones(N + 1)
    for p in primes:
        J[p::p] *= (p ** alpha - 1.0)
    return J


def divisor_sums(c, mu):
    """Y_d = sum_m c_{dm}, only for squarefree d (others are 0 since c lives on squarefree n)."""
    N = len(c) - 1
    Y = np.zeros(N + 1)
    for d in range(1, N + 1):
        if mu[d] != 0:
            Y[d] = c[d::d].sum()
    return Y


def main_terms(a, mu, primes, sigma):
    N = len(a) - 1
    n = np.arange(N + 1, dtype=float)
    n[0] = 1.0
    c1 = a / n ** (2 * sigma)
    c2 = a / n
    Y1 = divisor_sums(c1, mu)
    Y2 = divisor_sums(c2, mu)
    J1 = jordan_sqfree(N, 2 * sigma, primes)
    J2 = jordan_sqfree(N, 2 - 2 * sigma, primes)
    sq = mu != 0
    sq[0] = False
    M1 = float(mpmath.zeta(2 * sigma)) * float((J1[sq] * Y1[sq] ** 2).sum())
    S2 = float((J2[sq] * Y2[sq] ** 2).sum())
    return M1, S2


def m2_from_s2(S2, sigma, T):
    # zeta(2-2sigma) * (1/T) int_T^{2T} (t/2pi)^{1-2sigma} dt * S2
    e = 2 - 2 * sigma
    avg = (2 * np.pi / T) * (((2 * T) / (2 * np.pi)) ** e - (T / (2 * np.pi)) ** e) / e
    return float(mpmath.zeta(2 - 2 * sigma)) * avg * S2
