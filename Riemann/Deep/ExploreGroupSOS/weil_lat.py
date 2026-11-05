# weil_lat.py -- lattice (periodized) form of the truncated Weil functional, arbitrary period P >= 2L.
# Window: f in L^2(-L/2, L/2), L = log c (c real > 1). Q(f) = (1/P) sum_k m(2 pi k/P) |F(2 pi k/P)|^2,
# m(r) = W(e^{irx} 1_[-L,L]) = pole(r) + hplus(r) - T(r) - 2 sum_{n<c} w_n cos(r log n).
import math
import numpy as np
from scipy.special import psi

EULER = 0.5772156649015329

def vm(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return (p, math.log(p)) if m == 1 else None
    return None

def atoms(c, S=None):
    """prime powers 2 <= n < c (prime in S if S given): list of (log n, Lambda(n)/sqrt n)."""
    out = []
    n = 2
    while n < c:
        v = vm(n)
        if v is not None and (S is None or v[0] in S):
            out.append((math.log(n), v[1] / math.sqrt(n)))
        n += 1
    return out

def hplus(r):
    return psi(0.25 + 0.5j * np.asarray(r, dtype=float)).real - math.log(math.pi)

def Ttail(r, L, nterms=60):
    r = np.asarray(r, dtype=float)
    tot = np.zeros_like(r)
    for n in range(nterms):
        a = 2 * n + 0.5
        e = math.exp(-a * L)
        if e < 1e-300:
            break
        tot += e * (a * np.cos(r * L) - r * np.sin(r * L)) / (a * a + r * r)
    return -2 * tot

def pole(r, L):
    r = np.asarray(r, dtype=float)
    z1 = 1j * r + 0.5
    z2 = 1j * r - 0.5
    return (2 * np.sinh(z1 * L) / z1 + 2 * np.sinh(z2 * L) / z2).real

def prime(r, at):
    r = np.asarray(r, dtype=float)
    tot = np.zeros_like(r)
    for ln, w in at:
        tot += w * np.cos(r * ln)
    return -2 * tot

def m_sym(r, L, at, with_pole=True):
    v = hplus(r) - Ttail(r, L) + prime(r, at)
    if with_pole:
        v = v + pole(r, L)
    return v
