# Continuum limit as an ODE (zero pads, wb = R): r z' = z - T, T = (14 - 18r - 12z)/(30(3r-2)), z(7/9) = 0 (linear ODE).
# RK4 gain curve -> count balance -> infimum (golden in delta at x = 1/2, golden in b). Compare with PR141/143's exact contact.
import numpy as np, sys
from mpmath import mp, mpf, findroot
from system import Broot, Rd_native
from infimum import golden
def T(r, z): return (14 - 18 * r - 12 * z) / (30 * (3 * r - 2))
def f(r, z): return (z - T(r, z)) / r
def rk4(r0=7 / 9, r1=0.73, n=200000):
    h = (r1 - r0) / n; rs = [r0]; zs = [0.0]; r, z = r0, 0.0
    for _ in range(n):
        k1 = f(r, z); k2 = f(r + h / 2, z + h * k1 / 2); k3 = f(r + h / 2, z + h * k2 / 2); k4 = f(r + h, z + h * k3)
        z += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6; r += h; rs.append(r); zs.append(z)
    return np.array(rs[::-1]), np.array(zs[::-1])
rs, zs = rk4()
G = lambda w: np.interp(w, rs, zs)
def Rc(d, x):
    V = 5 / 6 - d; F = lambda c: (1 - d * (c + 2 * x * ((1 - c) / 2 + G(c)))) - (1 - d + V * (1.5 * c - 1))
    a, b = 0.735, 0.7775
    for _ in range(100):
        m = (a + b) / 2
        if (F(m) > 0) == (F(a) > 0): a = m
        else: b = m
    c = (a + b) / 2
    return min(1 - d + V * (1.5 * c - 1), Rd_native(d, x)), c
def Bb(b):   # max over delta at x = 1/2 of the root
    neg = lambda d: -Broot(b, d, 0.5, Rc(d, 0.5)[0])
    d = golden(neg, 0.38, 0.47, 1e-10)
    return -neg(d), d
if __name__ == '__main__':
    print("z(0.7524) = %.10f  z(0.7335) = %.10f" % (G(0.7524), G(0.7335)))
    bo = golden(lambda b: Bb(b)[0], 0.09, 0.13, 1e-11)
    B, d = Bb(bo); R, c = Rc(d, 0.5)
    print("ODE continuum contact: B* = %.13f  b = %.8f  delta = %.8f  crossing = %.8f  count = %.10f" % (B, bo, d, c, R))
    mp.dps = 30
    P = lambda t: 3888 * mpf(19)**5 * (15 * t - 13)**4 * (1 - t)**5 - 7 * (7 - 6 * t)**4 * (9 * t - 7)**5
    t = findroot(P, mpf('0.8733458'))
    print("PR143 exact root  = %s ; difference ODE - root = %.3e" % (mp.nstr(t, 22), B - float(t)))
