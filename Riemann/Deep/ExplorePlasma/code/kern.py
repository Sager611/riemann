# explore20: kernels for zeta (Phi), Davenport-Heilbronn (PhiDH) with closed-form derivatives
# Normalisation: Xi(z) = xi(1/2+iz) = 2*int_R Phi(u) e^{izu} du, Phi(u) = sum_n (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u}).
# DH: Lambda_f(1/2+iz) = 2*int_R PhiDH(u) e^{izu} du, PhiDH(u) = e^{3u/2} sum_n a_n n exp(-pi n^2 e^{2u}/5), a_n = 1, k, -k, -1, 0 (mod 5).
# Polymath dBN normalisation uses Phi_P(u) = Phi(2u); family Phi_P(u) e^{t u^2} <-> Phi(v) e^{s v^2}, t = 4 s.
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
sq5 = mp.sqrt(5)
kap = (mp.sqrt(10 - 2*sq5) - 2)/(sq5 - 1)
aDH = {0: 0, 1: 1, 2: kap, 3: -kap, 4: -1}
def terms_zeta(N=10):
    T = []
    for n in range(1, N+1):
        a = pi*n**2
        T += [(2*a**2, mp.mpf(9)/2, a), (-3*a, mp.mpf(5)/2, a)]
    return T
def terms_DH(N=30):
    return [(aDH[n % 5]*n, mp.mpf(3)/2, pi*n**2/5) for n in range(1, N+1) if n % 5]
TZ = terms_zeta(); TD = terms_DH()
def ev(T, u, k=0):
    """k-th derivative (k=0,1,2) of sum c exp(alpha u - a e^{2u}), u >= 0 (even kernels: use |u|)."""
    u = mp.mpf(u); x = mp.e**(2*u); s = mp.mpf(0)
    for c, al, a in T:
        e = c*mp.e**(al*u - a*x)
        if e == 0: continue
        if k == 0: s += e
        elif k == 1: s += (al - 2*a*x)*e
        else: s += ((al - 2*a*x)**2 - 4*a*x)*e
    return s
Phi = lambda u: ev(TZ, abs(u))
PhiDH = lambda u: ev(TD, abs(u))
def fDH(s):
    return 5**(-s)*mp.fsum(aDH[a]*mp.zeta(s, mp.mpf(a)/5) for a in range(1, 5))
def LamDH(s):
    return (5/pi)**((s+1)/2)*mp.gamma((s+1)/2)*fDH(s)
def xi(s):
    return mp.mpf(1)/2*s*(s-1)*pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
