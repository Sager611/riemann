# mp_core.py — Meixner–Pollaczek f_n = P_n^{(3/4)}(x; pi/2) (Romik 1902.06330), coefficients c_2n of Xi.
# f_n via (n+1) f_{n+1} = 2x f_n - (n+1/2) f_{n-1}, f_0 = 1, f_1 = 2x  (exact rationals).
from fractions import Fraction as Fr
import mpmath as mp

def f_coeffs(nmax):
    """list F[n] = monomial coefficient list (Fractions) of f_n(x), n <= nmax"""
    F = [[Fr(1)], [Fr(0), Fr(2)]]
    for n in range(1, nmax):
        a, b = F[n], F[n-1]
        c = [Fr(0)] * (n + 2)
        for k, v in enumerate(a):
            c[k+1] += 2 * v
        for k, v in enumerate(b):
            c[k] -= (n + Fr(1, 2)) * v
        F.append([v / (n + 1) for v in c])
    return F[:nmax+1]

def f_eval(nmax, x):
    """values f_0..f_nmax at x (mp number, real or complex) by the recurrence"""
    v = [mp.mpf(1), 2 * x]
    for n in range(1, nmax):
        v.append((2 * x * v[n] - (n + mp.mpf(1)/2) * v[n-1]) / (n + 1))
    return v[:nmax+1]

def f_hyp(n, x):
    """Romik's definition: (3/2)_n/n! i^n 2F1(-n, 3/4+ix; 3/2; 2)"""
    return mp.rf(mp.mpf(3)/2, n) / mp.factorial(n) * (1j)**n * mp.hyp2f1(-n, mp.mpf(3)/4 + 1j*x, mp.mpf(3)/2, 2)

def omega(x, M=12):
    s = mp.mpf(0)
    for n in range(1, M+1):
        s += (2*mp.pi**2*n**4*x**2 - 3*mp.pi*n**2*x) * mp.exp(-mp.pi*n**2*x)
    return s

def c2n_zeta(n):
    """c_2n = 4 sqrt2 int_1^inf omega(x) (x+1)^(-3/2) ((x-1)/(x+1))^(2n) dx, x = e^{2u}"""
    g = lambda u: omega(mp.exp(2*u)) * (mp.exp(2*u)+1)**(-mp.mpf(3)/2) * mp.tanh(u)**(2*n) * 2*mp.exp(2*u)
    return 4*mp.sqrt(2) * mp.quad(g, [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3.2])

def Xi(t):
    s = mp.mpf(1)/2 + 1j*t
    return (s*(s-1)/2 * mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)).real
