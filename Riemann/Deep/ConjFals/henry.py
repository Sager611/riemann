# Henry's conjecture (2601.00631 Conj; 2601.09276 Conj 1.1): for 0<x<1,
#   pi cot(pi x) + x < zeta(x) - psi(x) < pi cot(pi x) + (gamma+1/2) x + gamma - 1/2
from mpmath import mp, mpf, zeta, digamma, cot, pi, euler
mp.dps = 30
def DL(x): return zeta(x) - digamma(x) - pi*cot(pi*x) - x
def DU(x): return pi*cot(pi*x) + (euler+mpf(1)/2)*x + euler - mpf(1)/2 - (zeta(x) - digamma(x))
N=4000
mL=(1e9,None); mU=(1e9,None)
for k in range(1,N):
    x=mpf(k)/N
    a=DL(x); b=DU(x)
    if a<mL[0]: mL=(a,x)
    if b<mU[0]: mU=(b,x)
print('min D_L', mp.nstr(mL[0],12), 'at x=', mp.nstr(mL[1],8))
print('min D_U', mp.nstr(mU[0],12), 'at x=', mp.nstr(mU[1],8))
for x in (mpf('0.001'),mpf('0.01'),mpf('0.1'),mpf('0.5'),mpf('0.9'),mpf('0.99'),mpf('0.999')):
    print(mp.nstr(x,4), mp.nstr(DL(x),10), mp.nstr(DU(x),10))
