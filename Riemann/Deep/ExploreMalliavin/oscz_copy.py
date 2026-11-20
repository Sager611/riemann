# oscz.py (explore60): vectorized Hardy-type Z-functions on the critical line via Euler-Maclaurin Hurwitz zeta.
# Functions: 'zeta', 'dh' (Davenport-Heilbronn mod 5, odd, kappa as in ExplorePotential/dh_ft.py), 'chi8' (even, mod 8),
# 'chim8' (odd, mod 8). Z(t) = Re(e^{i theta(t)} L(1/2+it)); the discarded imaginary part is returned for checking.
import numpy as np
from scipy.special import loggamma
import mpmath as mp

PI = np.pi
B2 = np.array([float(mp.bernoulli(2*j)/mp.factorial(2*j)) for j in range(1, 13)])
KAP = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1)
# name -> (q, parity a, residue coefficients {r: c_r}); L(s) = q^{-s} sum_r c_r zeta(s, r/q)
FUN = {
    'zeta': (1, 0, {1: 1.0}),
    'dh': (5, 1, {1: 1.0, 2: KAP, 3: -KAP, 4: -1.0}),
    'chi8': (8, 0, {1: 1.0, 3: -1.0, 5: -1.0, 7: 1.0}),
    'chim8': (8, 1, {1: 1.0, 3: 1.0, 5: -1.0, 7: -1.0}),
}


def hurwitz(s, a, M, CH=16):  # explore62: CH 256 -> 16 (RSS grew with 10 MB temporaries)
    """zeta(s, a) for an array s, Euler-Maclaurin with M direct terms and 12 Bernoulli corrections."""
    kk = np.log(np.arange(M, dtype=float) + a)
    out = np.zeros(len(s), dtype=complex)
    for i in range(0, len(s), CH):
        ss = s[i:i + CH]
        out[i:i + CH] = np.exp(-ss[:, None]*kk[None, :]).sum(1)
    N = M + a
    out += N**(1 - s)/(s - 1) + 0.5*N**(-s)
    fac = s.copy()
    term = N**(-s - 1)
    for j in range(len(B2)):
        out += B2[j]*fac*term
        fac = fac*(s + 2*j + 1)*(s + 2*j + 2)
        term = term/(N*N)
    return out


def theta(name, t):
    q, a, _ = FUN[name]
    t = np.asarray(t, dtype=float)
    if name == 'zeta':
        return loggamma(0.25 + 0.5j*t).imag - 0.5*t*np.log(PI)
    return loggamma((0.5 + a)/2 + 0.5j*t).imag + 0.5*t*np.log(q/PI)


def Lval(name, s, M):
    q, a, co = FUN[name]
    out = np.zeros(len(s), dtype=complex)
    for r, c in co.items():
        out += c*hurwitz(s, r/q, M)
    return out*q**(-s) if q > 1 else out


def Z(name, t, chunk=4000):
    """Returns (Z, imag part of e^{i theta} L) on the array t (needs t >= 0)."""
    t = np.asarray(t, dtype=float)
    z = np.empty(len(t))
    im = np.empty(len(t))
    for i in range(0, len(t), chunk):
        tt = t[i:i + chunk]
        M = int(tt.max()/2) + 60
        v = np.exp(1j*theta(name, tt))*Lval(name, 0.5 + 1j*tt, M)
        z[i:i + chunk] = v.real
        im[i:i + chunk] = v.imag
    return z, im


def Z_mp(name, t):
    """mpmath reference value (for spot checks)."""
    q, a, co = FUN[name]
    s = mp.mpf(0.5) + 1j*mp.mpf(t)
    if name == 'zeta':
        return float(mp.siegelz(t))
    L = mp.power(q, -s)*sum(c*mp.zeta(s, mp.mpf(r)/q) for r, c in co.items())
    return float(mp.re(mp.exp(1j*mp.mpf(theta(name, np.array([t]))[0]))*L))
