"""Lattice symbol mu(t) of Weil's form on [-L/2, L/2], L = ln c, at ANY real t (float64).
mu(t) = phi(t) + 2 sum_n c^{-b_n}[b_n cos(tL) - t sin(tL)]/(b_n^2+t^2) - 2 sum_{m<c} Lambda(m) m^{-1/2} cos(t log m),
phi(t) = -log pi + Re psi(1/4 + i t/2), b_n = 2n + 1/2.  At t = pi k/L this is Chekmachev's mu(k) (sin(pi k) = 0)."""
import math
import numpy as np

def vonmangoldt_table(c):
    out = []
    for m in range(2, int(math.ceil(c))):
        if m >= c:
            break
        n, p = m, None
        for q in range(2, m + 1):
            if n % q == 0:
                p = q
                break
        while n % p == 0:
            n //= p
        if n == 1:
            out.append((m, math.log(p)))
    return out

def re_digamma_quarter(t):
    """Re psi(1/4 + i t/2), float64, vectorized; recurrence up to |z| >= 12 then asymptotic series."""
    t = np.asarray(t, dtype=float)
    acc = np.zeros_like(t)
    x = np.full_like(t, 0.25)
    y = t / 2.0
    # psi(z) = psi(z+N) - sum_{j<N} 1/(z+j); use N = 12 when |y| small
    small = np.abs(y) < 12.0
    N = 12
    if np.any(small):
        for j in range(N):
            xr = x[small] + j
            acc[small] -= xr / (xr * xr + y[small] ** 2)
        x = x.copy()
        x[small] += N
    z = x + 1j * y
    zi = 1.0 / z
    z2 = zi * zi
    ser = np.log(z) - 0.5 * zi - z2 * (1.0 / 12 - z2 * (1.0 / 120 - z2 * (1.0 / 252 - z2 * (1.0 / 240 - z2 * (1.0 / 132)))))
    return acc + ser.real

def mu(t, c, nb=60):
    t = np.asarray(t, dtype=float)
    L = math.log(c)
    val = -math.log(math.pi) + re_digamma_quarter(t)
    cl, sl = np.cos(t * L), np.sin(t * L)
    for n in range(nb):
        b = 2 * n + 0.5
        w = c ** (-b)
        if w < 1e-300:
            break
        val += 2 * w * (b * cl - t * sl) / (b * b + t * t)
    for m, lam in vonmangoldt_table(c):
        val -= 2 * lam / math.sqrt(m) * np.cos(t * math.log(m))
    return val

def sigma(c):
    return sum(l / math.sqrt(m) for m, l in vonmangoldt_table(c))
