
import mpmath as mp, sys
mp.mp.dps = 25
pi = mp.pi
def phi(u):
    e2 = mp.e**(2*u)
    return sum((2*pi**2*k**4*mp.e**(mp.mpf(9)/2*u)-3*pi*k**2*mp.e**(mp.mpf(5)/2*u))*mp.e**(-pi*k**2*e2) for k in range(1,5))
def ustar(n):
    return mp.findroot(lambda u: u*(pi*mp.e**(2*u)-mp.mpf(9)/4)-n, mp.log(n)/2-1)
def Mt(k, lo=0, hi=5, pts=11):
    return mp.quad(lambda u: u**(2*k)*phi(u), mp.linspace(lo, hi, pts))
n = 3
gl = {}
def gam(m):
    if m not in gl: gl[m] = 32*mp.factorial(m)*Mt(m)/mp.factorial(2*m)
    return gl[m]
w = mp.mpc(-50, 10)
direct = sum(gam(n+k)*w**k/mp.factorial(k) for k in range(0, 60))
mix = gam(n)*mp.quad(lambda u: u**(2*n)*phi(u)*mp.hyp0f1(n+mp.mpf(1)/2, u**2*w/4), mp.linspace(0,5,11))/Mt(n)
print('identity n=3 w=-50+10i: direct', mp.nstr(direct,12), ' mixture', mp.nstr(mix,12), flush=True)
