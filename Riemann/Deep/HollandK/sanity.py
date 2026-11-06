#!/usr/bin/env python3
# Sanity checks (not part of the proof): normalization of the kernel, the saddle main term vs quadrature,
# Jensen hyperbolicity for small d, n (Sturm), and the parameter pipeline at n ~ 1e45 via the main term.
import mpmath as mp, sys
mp.mp.dps = 60
pi = mp.pi
def Phi(u):   # Holland's Phi = 2 * Polya kernel: Xi(t) = int_0^oo Phi(u) cos(t u) du
    s = mp.mpf(0)
    for k in range(1, 8):
        s += (2*pi**2*k**4*mp.exp(9*u) - 3*pi*k**2*mp.exp(5*u))*mp.exp(-pi*k**2*mp.exp(4*u))
    return 2*s
# 1) Xi(0) = xi(1/2)
xi_half = mp.mpf('0.497120778188314109912')   # xi(1/2) = 0.4971207781883141...
I0 = mp.quad(Phi, [0, 0.2, 0.5, 1, 2])
print("1) int Phi =", mp.nstr(I0, 20), " xi(1/2) =", mp.nstr(xi_half, 20))
# gamma(n) = n!/(2n)! M_n, check against Taylor coefficients of xi(1/2+z) via mpmath differentiation for n = 0,1,2
def M(n):
    f = lambda u: Phi(u)*u**(2*n)
    return mp.quad(f, [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2])
xi = lambda s: s*(s-1)/2*pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
for n in range(3):
    g_quad = mp.factorial(n)/mp.factorial(2*n)*M(n)
    g_tay = mp.diff(lambda z: xi(mp.mpf(1)/2+z), 0, 2*n)/mp.factorial(2*n)*mp.factorial(n)
    print("   gamma(%d): moment %s  Taylor %s" % (n, mp.nstr(g_quad, 15), mp.nstr(g_tay, 15)))
# 2) saddle main term for Mt(z) = int_0^oo phi(u) u^{2z} du, phi(u) = sum (2 pi^2 k^4 e^{9u/2} - 3 pi k^2 e^{5u/2}) e^{-pi k^2 e^{2u}}
def u0(z):
    u = mp.log(z)/2
    for _ in range(200):
        u = mp.log((z/u + mp.mpf(9)/4)/pi)/2
    return u
def main(z):
    u = u0(z)
    f = 2*z*mp.log(u) + 9*u/2 - pi*mp.exp(2*u)
    fpp = -2*z/u**2 - 4*pi*mp.exp(2*u)
    return 2*pi**2*mp.exp(f)*mp.sqrt(2*pi/(-fpp))
def Mt(z, U):
    def g(v):
        s = 0
        for k in range(1, 4):
            s += (2*pi**2*k**4*mp.exp(9*v/2) - 3*pi*k**2*mp.exp(5*v/2))*mp.exp(-pi*k**2*mp.exp(2*v))
        return s*mp.exp(2*z*mp.log(v))
    return mp.quad(g, [U*0.5, U*0.9, U-0.3, U-0.1, U, U+0.1, U+0.3, U+0.6, U+1.5])
for x0 in (200, 2000):
    for z in (mp.mpf(x0), mp.mpc(x0, x0/10), mp.mpc(x0*0.95, x0/12)):
        U = mp.re(u0(z))
        q = Mt(z, U)/main(z)
        print("2) z =", mp.nstr(z, 8), " Mt/main - 1 =", mp.nstr(q-1, 6), " e^{-Re u0} =", mp.nstr(mp.exp(-U), 4))
