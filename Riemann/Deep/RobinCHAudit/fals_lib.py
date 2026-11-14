#!/usr/bin/env python3
# fals.py (robinchaudit) -- falsification tests of RobinCH Lemma 11_k ((B2_k), (B1'), m_c = E|V - c|, MAD_k) and Lemma 13
# at corner parameters, by a method different from RobinCH/test_b2k.py: F_phi(rho) = E F(rho; L' + V) expanded in Taylor
# series of G(y) = e^{-sy} F(rho; y) around L' (exact E_n closed forms) against the moments E[V^n e^{sV}] = d^n/ds^n M(s),
# M(s) = ((e^{sa} - 1)/(sa))^k. Non-rigorous (mpmath, 60 digits); a ratio true/bound > 1 would falsify.
#   uv run --no-project --with mpmath python3 fals.py
import itertools, math, time
import mpmath as mp
mp.mp.dps = 60


def Fy(s, y):       # F(rho; y) = int_y^inf e^{su} (1/u + 1/u^2) du = E_1(-sy) + E_2(-sy)/y
    z = -s*y
    return mp.expint(1, z) + mp.expint(2, z)/y

def Ij(j, z):        # m^{(j)}(z) = int_0^1 t^j e^{zt} dt
    if abs(z) < 1:
        return sum(z**i/(mp.factorial(i)*(i + j + 1)) for i in range(90))
    v = (mp.exp(z) - 1)/z
    for i in range(1, j + 1): v = (mp.exp(z) - i*v)/z
    return v

def Mders(s, a, k, N):   # [M^{(n)}(s) for n <= N], M(s) = m(sa)^k: Taylor series of m around z0 = sa, power k, exact
    z0 = s*a
    c = [Ij(j, z0)/mp.factorial(j) for j in range(N + 1)]
    p = [mp.mpc(1)] + [mp.mpc(0)]*N
    for _ in range(k):
        p = [sum(p[i]*c[n - i] for i in range(n + 1)) for n in range(N + 1)]
    return [p[n]*mp.factorial(n)*a**n for n in range(N + 1)]

def Fphi(s, Lp, lam, k, N=10):   # E F(rho; L'+V) = sum_n (-1)^n [L'^-n E_{n+1}(-sL') + (n+1) L'^{-n-1} E_{n+2}(-sL')] M^{(n)}(s)
    a = mp.mpf(lam)/k
    Md = Mders(s, a, k, N)
    z = -s*Lp; tot = mp.mpc(0); last = 0
    for n in range(N + 1):
        t = (-1)**n*(Lp**(-n)*mp.expint(n + 1, z) + (n + 1)*Lp**(-n - 1)*mp.expint(n + 2, z))*Md[n]/mp.factorial(n)
        tot += t; last = t
    return tot, last

def MADk(k):
    return 2*sum((-1)**j*math.comb(k, j)*(mp.mpf(k)/2 - j)**(k + 1) for j in range(k//2 + 1))/mp.factorial(k + 1)

