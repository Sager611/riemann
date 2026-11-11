
import mpmath as mp
mp.mp.dps = 30
def xi(s): return 0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
co = mp.taylor(lambda z: 8*xi(mp.mpf(1)/2+z), 0, 12)
gam = [co[2*n]*mp.factorial(n) for n in range(7)]
def phi(u):
    return mp.nsum(lambda k: (2*mp.pi**2*k**4*mp.e**(4.5*u)-3*mp.pi*k**2*mp.e**(2.5*u))*mp.e**(-mp.pi*k**2*mp.e**(2*u)), [1, mp.inf])
def Mt(k): return mp.quad(lambda u: u**(2*k)*phi(u), [0, 0.5, 1, 2, 4])
for n in range(7):
    m = Mt(n)
    print(n, mp.nstr(gam[n],15), mp.nstr(gam[n]/(mp.factorial(n)*m/mp.factorial(2*n)),15))
print('odd coeffs', [mp.nstr(co[2*n+1],3) for n in range(3)])
