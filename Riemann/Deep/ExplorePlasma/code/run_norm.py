import sys, time; sys.path.insert(0, __file__.rsplit('/', 1)[0])
# normalisation checks: Xi = 2 int Phi e^{izu} = xi(1/2+iz); DH kernel evenness (functional equation) and transform; DH off-line zero
from kern import *
mp.mp.dps = 30
def PhiDH_raw(u, N=80):  # no |u|: tests evenness = functional equation with root number 1
    u = mp.mpf(u); x = mp.e**(2*u)
    return mp.e**(3*u/2)*mp.fsum(aDH[n % 5]*n*mp.e**(-pi*n**2*x/5) for n in range(1, N+1))
print('Xi(0): kernel', mp.nstr(2*mp.quad(Phi, [-3, -1, 0, 1, 3]), 15), ' xi(1/2)', mp.nstr(xi(mp.mpf(0.5)), 15))
for z in [5, 10]:
    print('Xi(%d): kernel' % z, mp.nstr(4*mp.quad(lambda u: Phi(u)*mp.cos(z*u), [0, 0.5, 1, 2, 3]), 12), ' xi', mp.nstr(mp.re(xi(mp.mpf(0.5) + 1j*z)), 12))
print('kappa =', mp.nstr(kap, 15))
for u in [0.05, 0.2, 0.4]:
    print('PhiDH(+%.2f) = %s   PhiDH(-%.2f) = %s' % (u, mp.nstr(PhiDH_raw(u), 15), u, mp.nstr(PhiDH_raw(-u), 15)))
print('XiDH(0): kernel', mp.nstr(2*mp.quad(PhiDH, [-3, -1, 0, 1, 3]), 15), ' Lambda_DH(1/2)', mp.nstr(LamDH(mp.mpf(0.5)), 15))
print('XiDH(10): kernel', mp.nstr(4*mp.quad(lambda u: PhiDH(u)*mp.cos(10*u), [0, 0.5, 1, 2, 3]), 12), ' Lambda_DH(1/2+10i)', mp.nstr(mp.re(LamDH(0.5 + 10j)), 12))
print('FE: Lambda_DH(0.3+2i) =', mp.nstr(LamDH(0.3 + 2j), 10), '  conj Lambda_DH(0.7+2i) =', mp.nstr(mp.conj(LamDH(0.7 + 2j)), 10))
s0 = mp.mpc('0.808517', '85.699348')
print('DH off-line zero s0 = 0.808517+85.699348i: |Lambda_DH(s0)| =', mp.nstr(abs(LamDH(s0)), 4), '  |Lambda_DH(1/2+85.699348i)| =', mp.nstr(abs(LamDH(0.5 + 85.699348j)), 4))
