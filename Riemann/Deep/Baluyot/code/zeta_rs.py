"""Vectorized zeta(sigma+it) for large t: Riemann-Siegel main sums + leading correction.

zeta(s) = sum_{n<=m} n^{-s} + chi(s) sum_{n<=m} n^{s-1} + R,  m = floor(sqrt(t/2pi)),
R ~ (-1)^{m-1} (t/2pi)^{-sigma/2} e^{-i theta(t)} Psi(tau),  tau = sqrt(t/2pi) - m,
Psi(tau) = cos(2pi(tau^2 - tau - 1/16)) / cos(2pi tau).
Accuracy is measured against python-flint (Arb) in check_rs.py, not assumed.
"""
import numpy as np
from scipy.special import loggamma

LOG2PI = np.log(2 * np.pi)


def theta(t):
    t = np.asarray(t, dtype=float)
    return np.imag(loggamma(0.25 + 0.5j * t)) - 0.5 * t * np.log(np.pi)


def log_chi(sigma, t):
    # chi(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s); log via loggamma, sin via exp form
    s = sigma + 1j * np.asarray(t, dtype=float)
    # log sin(pi s/2) for large Im s: sin(z) = (e^{iz}-e^{-iz})/(2i), Im z = pi t/2 > 0 dominant e^{-iz}
    z = np.pi * s / 2
    logsin = -1j * z + np.log((np.exp(2j * z) - 1) / 2j)
    return s * np.log(2) + (s - 1) * np.log(np.pi) + logsin + loggamma(1 - s)


def psi(tau):
    return np.cos(2 * np.pi * (tau * tau - tau - 1.0 / 16)) / np.cos(2 * np.pi * tau)


def zeta_rs(sigma, t):
    t = np.asarray(t, dtype=float)
    a = np.sqrt(t / (2 * np.pi))
    m = np.floor(a).astype(np.int64)
    mmax = int(m.max())
    n = np.arange(1, mmax + 1, dtype=float)
    logn = np.log(n)
    mask = n[None, :] <= m[:, None]
    ph = np.exp(-1j * np.outer(t, logn))  # n^{-it}
    s1 = (ph * (n ** (-sigma))[None, :] * mask).sum(axis=1)
    s2 = (np.conj(ph) * (n ** (sigma - 1))[None, :] * mask).sum(axis=1)
    chi = np.exp(log_chi(sigma, t))
    tau = a - m
    R = ((-1.0) ** (m - 1)) * (t / (2 * np.pi)) ** (-sigma / 2) * np.exp(-1j * theta(t)) * psi(tau)
    return s1 + chi * s2 + R


def zeta_rs_blocks(sigma, t, block=512):
    out = np.empty(len(t), dtype=complex)
    for i in range(0, len(t), block):
        out[i:i + block] = zeta_rs(sigma, t[i:i + block])
    return out
