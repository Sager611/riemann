#!/usr/bin/env python3
"""Numerical checks for the angular (k-th derivative) theta transformation used in NOTE.md.
1) Mellin-Bessel: int_0^inf v^(mu-1) K_{1/3}(a v) dv = 2^(mu-2) a^(-mu) G((mu+1/3)/2) G((mu-1/3)/2), mu = 2s+k.
2) Chain rule at the cusp map (A.9) of OpenAI 5 Oct 2026:
   z' = -d/c - zbar/(c^2 (v^2+|z|^2)),  v' = v/(|c|^2 (v^2+|z|^2));
   claim: d_zbar^k [Phi(z',v')] at z=0 equals (-1/(c^2 v^2))^k (d_w^k Phi)(z'_0, v'_0), here k = 2, 3,
   for Phi(w,u) = e(l w) u K_{1/3}(4 pi |l| u), e(l w) = exp(2 pi i (l w + conj(l w))).
Run: python3 angular_check.py   (mpmath only; ~20 s, <40 MB)
"""
import mpmath as mp

mp.mp.dps = 30
print("1) Mellin-Bessel identity, relative errors:")
for s, k in ((0.7, 2), (1.3, 3), (0.9+0.4j, 4)):
    mu = 2*s+k
    a = 4*mp.pi*1.7
    lhs = mp.quad(lambda v: v**(mu-1)*mp.besselk(mp.mpf(1)/3, a*v), [0, 1, mp.inf])
    rhs = 2**(mu-2)*a**(-mu)*mp.gamma((mu+mp.mpf(1)/3)/2)*mp.gamma((mu-mp.mpf(1)/3)/2)
    print(f"   s={s}, k={k}: {mp.nstr(abs(lhs/rhs-1), 3)}")

c = mp.mpc(2, 1); d = mp.mpc(0.3, -0.7); v0 = mp.mpf("0.6"); l = mp.mpc(0.5, 1.2)

def phi(w, u):
    return mp.exp(2j*mp.pi*(l*w+mp.conj(l*w)))*u*mp.besselk(mp.mpf(1)/3, 4*mp.pi*abs(l)*u)

def F(x, y):
    z = mp.mpc(x, y); r2 = v0**2+abs(z)**2
    return phi(-d/c-mp.conj(z)/(c**2*r2), v0/(abs(c)**2*r2))

print("2) chain rule, relative errors:")
w0, u0 = -d/c, 1/(abs(c)**2*v0)
for k in (2, 3):
    # d_zbar^k = 2^-k (d_x + i d_y)^k, expanded binomially into mixed partials
    tot = 0
    for j in range(k+1):
        tot += mp.binomial(k, j)*(1j)**j*mp.diff(F, (0, 0), (k-j, j))
    lhs = tot/2**k
    rhs = (-1/(c**2*v0**2))**k*(2j*mp.pi*l)**k*phi(w0, u0)
    print(f"   k={k}: {mp.nstr(abs(lhs/rhs-1), 3)}")
