#!/usr/bin/env python3
"""Numerical checks for the angular (k-th derivative) theta transformation used in NOTE.md (k = 1 + m >= 1).
Notation of OpenAI, 5 Oct 2026 (paper2), App. A.2: e(lw) = exp(2 pi i (l w + conj(l w))), lambda = sqrt(-3),
(A.9): g^{-1}(z + a/c, v) = (-d/c - zbar/(c^2 (v^2+|z|^2)), v/(|c|^2 (v^2+|z|^2))).
1) Mellin-Bessel: int_0^inf v^(mu-1) K_{1/3}(a v) dv = 2^(mu-2) a^(-mu) G((mu+1/3)/2) G((mu-1/3)/2), mu = 2s+k.
2) Lemma A (chain rule): d_zbar^k [Phi(g^{-1}(z+a/c,v))] at z=0 = (-1/(c^2 v^2))^k (d_w^k Phi)(-d/c, 1/(|c|^2 v)),
   k = 1..4, for a Fourier mode Phi(w,u) = e(l w) u K_{1/3}(4 pi |l| u) and for the constant term Phi0 = u^(2/3)
   (whose right side is 0: the constant terms are killed).
3) Mellin constants for one mode: forward  int_0^inf (2 pi i lbar)^k v K(4pi|l|v) v^(2s+k-2) dv
   = (i^k/4)(2pi)^(-2s) |l|^(-2s) (lbar/|l|)^k G(s+k/2+1/6) G(s+k/2-1/6); with l = lambda^-3 m this is
   (1/4)(27/(2pi)^2)^s N(m)^-s conj(m/|m|)^k G G  [c_k = 3^(5/2)/4 for every k];  dual (d_z^k, l instead of lbar):
   (i^k/4)(2pi)^(-2s) N(l)^-s (l/|l|)^k G G   [the dual factor is (2pi)^(-2s), not (27/(2pi)^2)^s].
4) Functional equation of one translated mode (v -> 1/(|c|^2 u)):
   int_0^inf d_zbar^k[Phi(g^{-1}(z+a/c,v))]_{z=0} v^(2s+k-2) dv
   = (-1)^k (cbar/c)^k |c|^(2-4s) (i^k/4)(2pi)^(-2s') N(l)^-s' (l/|l|)^k e(-d l/c) G(s'+k/2+1/6) G(s'+k/2-1/6), s'=1-s.
Run: python3 angular_check.py   (mpmath only; measured 16 s, 22 MB peak RSS)
"""
import mpmath as mp

mp.mp.dps = 30
third = mp.mpf(1)/3
G = mp.gamma
print("1) Mellin-Bessel identity, relative errors:")
for s, k in ((0.7, 2), (1.3, 3), (0.9+0.4j, 4)):
    mu = 2*s+k
    a = 4*mp.pi*1.7
    lhs = mp.quad(lambda v: v**(mu-1)*mp.besselk(third, a*v), [0, 1, mp.inf])
    rhs = 2**(mu-2)*a**(-mu)*G((mu+third)/2)*G((mu-third)/2)
    print(f"   s={s}, k={k}: {mp.nstr(abs(lhs/rhs-1), 3)}")

c = mp.mpc(2, 1); d = mp.mpc(0.3, -0.7); v0 = mp.mpf("0.6"); l = mp.mpc(0.5, 1.2)
N = lambda x: abs(x)**2
ebr = lambda x: mp.exp(2j*mp.pi*(x+mp.conj(x)))


def phi(w, u):
    return ebr(l*w)*u*mp.besselk(third, 4*mp.pi*abs(l)*u)


def phi0(w, u):
    return u**(mp.mpf(2)/3)


def dzbar_k(Phi, k):
    """d_zbar^k of z -> Phi(g^{-1}(z+a/c, v0)) at z = 0, via d_zbar = (d_x + i d_y)/2 (binomial mixed partials)."""
    def F(x, y):
        z = mp.mpc(x, y); r2 = v0**2+abs(z)**2
        return Phi(-d/c-mp.conj(z)/(c**2*r2), v0/(N(c)*r2))
    tot = 0
    for j in range(k+1):
        tot += mp.binomial(k, j)*(1j)**j*mp.diff(F, (0, 0), (k-j, j))
    return tot/2**k


print("2) Lemma A: relative error for the mode, absolute value for the constant term (should be ~0):")
w0, u0 = -d/c, 1/(N(c)*v0)
for k in (1, 2, 3, 4):
    lhs = dzbar_k(phi, k)
    rhs = (-1/(c**2*v0**2))**k*(2j*mp.pi*l)**k*phi(w0, u0)
    cst = dzbar_k(phi0, k)
    print(f"   k={k}: mode {mp.nstr(abs(lhs/rhs-1), 3)}   constant term {mp.nstr(abs(cst), 3)}")

print("3) one-mode Mellin constants, relative errors (forward generic l, forward l = lambda^-3 m, dual):")
lam = mp.sqrt(-3)
m = mp.mpc(4, 3*mp.sqrt(3)/1)  # any nonzero m; only |m| and arg m enter
for s, k in ((mp.mpc(1.3, 0.4), 1), (mp.mpc(1.1, -0.7), 2), (mp.mpc(1.6, 0.2), 3)):
    gg = G(s+k/mp.mpf(2)+mp.mpf(1)/6)*G(s+k/mp.mpf(2)-mp.mpf(1)/6)
    def fwd(ll):
        return mp.quad(lambda v: (2j*mp.pi*mp.conj(ll))**k*v*mp.besselk(third, 4*mp.pi*abs(ll)*v)*v**(2*s+k-2), [0, 1, mp.inf])
    r1 = (1j**k/4)*(2*mp.pi)**(-2*s)*abs(l)**(-2*s)*(mp.conj(l)/abs(l))**k*gg
    ll = m/lam**3
    r2 = mp.mpf(1)/4*(27/(2*mp.pi)**2)**s*N(m)**(-s)*(mp.conj(m)/abs(m))**k*gg
    dual = mp.quad(lambda v: (2j*mp.pi*l)**k*v*mp.besselk(third, 4*mp.pi*abs(l)*v)*v**(2*s+k-2), [0, 1, mp.inf])
    r3 = (1j**k/4)*(2*mp.pi)**(-2*s)*N(l)**(-s)*(l/abs(l))**k*gg
    print(f"   k={k}: {mp.nstr(abs(fwd(l)/r1-1), 3)}  {mp.nstr(abs(fwd(ll)/r2-1), 3)}  {mp.nstr(abs(dual/r3-1), 3)}")

print("4) functional equation of one translated mode, relative error:")
for s, k in ((mp.mpc(0.3, 0.2), 1), (mp.mpc(-0.4, 0.5), 2), (mp.mpc(0.1, -0.3), 3)):
    sp = 1-s
    def integrand(v):
        u = 1/(N(c)*v)
        dw = (2j*mp.pi*l)**k*ebr(-d*l/c)*u*mp.besselk(third, 4*mp.pi*abs(l)*u)  # (d_w^k Phi)(-d/c, u)
        return (-1/(c**2*v**2))**k*dw*v**(2*s+k-2)
    lhs = mp.quad(integrand, [0, mp.mpf(1)/N(c), 1, mp.inf])
    gg = G(sp+k/mp.mpf(2)+mp.mpf(1)/6)*G(sp+k/mp.mpf(2)-mp.mpf(1)/6)
    rhs = (-1)**k*(mp.conj(c)/c)**k*N(c)**(1-2*s)*(1j**k/4)*(2*mp.pi)**(-2*sp)*N(l)**(-sp)*(l/abs(l))**k*ebr(-d*l/c)*gg
    print(f"   s={mp.nstr(s, 3)}, k={k}: {mp.nstr(abs(lhs/rhs-1), 3)}")
