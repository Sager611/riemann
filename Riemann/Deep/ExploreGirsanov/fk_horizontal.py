# fk_horizontal.py (explore56): Feynman-Kac (Gaussian) form of the de Bruijn-Newman flow.
# H_t(z) = int_0^inf e^{t u^2} Phi(u) cos(z u) du (Polymath15 normalisation, H_0(z) = xi(1/2 + i z/2)/8).
# Claim (t >= 0): H_t(z) = (1/8) E[ xi(1/2 + i z/2 + sqrt(t/2) G) ], G ~ N(0,1): the heat semigroup in sigma = Re s.
# (Hubbard-Stratonovich e^{t u^2} = E e^{sqrt(2t) G u} and xi(1/2 + w) = 8 int Phi(u) cosh(2 w u) du.)
import mpmath as mp
mp.mp.dps = 25
def Phi(u, N=12):
    return sum((2*mp.pi**2*n**4*mp.e**(9*u) - 3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u)) for n in range(1, N+1))
def H_kernel(t, z):
    return mp.quad(lambda u: mp.e**(t*u*u)*Phi(u)*mp.cos(z*u), [0, 0.5, 1, 2, 3])
def xi(s):
    return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def H_fk(t, z):
    if t == 0: return xi(mp.mpf(1)/2 + 1j*z/2)/8
    c = mp.sqrt(mp.mpf(t)/2)
    f = lambda g: xi(mp.mpf(1)/2 + 1j*mp.mpf(z)/2 + c*g)*mp.e**(-g*g/2)/mp.sqrt(2*mp.pi)
    return mp.quad(f, [-14, -6, -2, 0, 2, 6, 14])/8
for t in [0, 0.05, 0.1, 0.2]:
    for z in [0, 10, 28.3, 40, 60]:
        a = H_kernel(mp.mpf(t), mp.mpf(z)); b = H_fk(t, z)
        print('t=%-5g z=%-5g kernel %s  FK %s  |diff| %s' % (t, z, mp.nstr(a, 12), mp.nstr(b, 12), mp.nstr(abs(a-b), 3)))
# Gaussian peak shift: weight |xi(sigma + iT)| exp(-(sigma-mu)^2/t) at T = 1e6 peaks near mu + t*L/4, L = log(T/2pi)
T = mp.mpf(10)**6; L = mp.log(T/(2*mp.pi)); t = mp.mpf('0.1'); mu = mp.mpf('0.4')
sig = [mu + mp.mpf(k)/50 for k in range(-25, 101)]
w = [abs(xi(s + 1j*T))*mp.e**(-(s - mu)**2/t) for s in sig]
k = max(range(len(w)), key=lambda i: w[i])
print('T=1e6, t=0.1, mu=0.4: weight peak at sigma = %s, mu + tL/4 = %s' % (mp.nstr(sig[k], 4), mp.nstr(mu + t*L/4, 4)))
