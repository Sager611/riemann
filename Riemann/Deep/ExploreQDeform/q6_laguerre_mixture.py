#!/usr/bin/env python3
"""q6_laguerre_mixture.py -- explore16: numeric check of the side identity
J^{d,n}(X) = 8 d!/(4^n (1/2)_n (n+1/2)_d) * int_R Phi(u) u^{2n} L_d^{(n-1/2)}(-u^2 X/4) du,
Phi normalised by xi(1/2 + z) = int_R Phi(u) cosh(z u) du, i.e. Phi(u) = 2 sum_n (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u}).
Run: uv run --no-project --with mpmath python3 -u q6_laguerre_mixture.py
"""
import mpmath as mp
mp.mp.dps = 40

def Phi(u):
    s = mp.mpf(0)
    for n in range(1, 8):
        s += 2 * (2 * mp.pi**2 * n**4 * mp.e**(9 * u / 2) - 3 * mp.pi * n**2 * mp.e**(5 * u / 2)) * mp.e**(-mp.pi * n**2 * mp.e**(2 * u))
    return s

def b(m):
    return 2 * mp.quad(lambda u: Phi(u) * u**(2 * m), [0, 0.5, 1, 2, 4])

def main():
    gam = [8 * mp.factorial(m) * b(m) / mp.factorial(2 * m) for m in range(0, 9)]
    xi_half = mp.mpf(1) / 2 * mp.mpf(1) / 2 * (-mp.mpf(1) / 2) * mp.pi**(-mp.mpf(1) / 4) * mp.gamma(mp.mpf(1) / 4) * mp.zeta(mp.mpf(1) / 2)
    print('gamma(0)/8 =', mp.nstr(gam[0] / 8, 15), ' xi(1/2) =', mp.nstr(xi_half, 15))
    for (d, n, X) in ((3, 2, mp.mpf('-0.7')), (4, 1, mp.mpf('-2.5')), (5, 3, mp.mpf('-11.0'))):
        J = sum(mp.binomial(d, j) * gam[n + j] * X**j for j in range(d + 1))
        a = n - mp.mpf(1) / 2
        pref = 8 * mp.factorial(d) / (4**n * mp.rf(mp.mpf(1) / 2, n) * mp.rf(n + mp.mpf(1) / 2, d))
        M = pref * 2 * mp.quad(lambda u: Phi(u) * u**(2 * n) * mp.laguerre(d, a, -u**2 * X / 4), [0, 0.5, 1, 2, 4])
        print(f'd={d} n={n} X={mp.nstr(X, 3)}: J = {mp.nstr(J, 15)}, mixture = {mp.nstr(M, 15)}, rel diff = {mp.nstr(abs(J - M) / abs(J), 3)}')

if __name__ == '__main__':
    main()
