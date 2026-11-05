# t2_sensitivity.py — strange-identity coefficients vs zeros: (a) FE link L(-2k-1,chi) <-> L(2k+2,chi) (Glaisher numbers);
# (b) Stokes constants = Dirichlet coefficients, peeled from the exact sequence; (c) relative change of the coefficients if the
# lowest zero of L(s,chi12) were moved off the line by delta (FE-preserving quadruple), vs. integrality.
import mpmath as mp
from fractions import Fraction as Fr
from math import comb
mp.mp.dps = 400
def chi12(n): return {1: 1, 11: 1, 5: -1, 7: -1}.get(n % 12, 0)
def bern(n):
    B = [Fr(0)]*(n+1); B[0] = Fr(1)
    for m in range(1, n+1): B[m] = -sum(comb(m+1, k)*B[k] for k in range(m))/(m+1)
    return B
BN = bern(300)
def Bpoly(n, x): return sum(comb(n, k)*BN[k]*x**(n-k) for k in range(n+1))
def Lneg(m, chi, f): return -Fr(f)**m/(m+1)*sum(chi(r)*Bpoly(m+1, Fr(r, f)) for r in range(1, f+1))
f = 12
print('Glaisher T_k = -L(-2k-1,chi12)/2 * (-1)^k:', [int(Lneg(2*k+1, chi12, f)*(-1)**(k+1)/2) for k in range(8)])
# (a) FE: L(-2k-1,chi) = (f/pi)^(2k+3/2) * k! * Gamma(k+3/2) * (-1)^(k+1) / pi * L(2k+2,chi)   (chi real even primitive, eps=1)
L = lambda s: mp.dirichlet(s, [chi12(n) for n in range(12)])
for k in [0, 1, 5, 20]:
    lhs = mp.mpf(Lneg(2*k+1, chi12, f).numerator)/Lneg(2*k+1, chi12, f).denominator
    rhs = (mp.mpf(f)/mp.pi)**(2*k+mp.mpf(3)/2)*mp.factorial(k)*mp.gamma(k+mp.mpf(3)/2)*(-1)**(k+1)/mp.pi*L(2*k+2)
    print(f'k={k}: L(-2k-1) exact / FE-from-L(2k+2) = {mp.nstr(lhs/rhs, 20)}')
# (b) peel Stokes constants from the exact coefficient sequence alone
def Lpos(k):  # L(2k+2,chi) recovered from the exact coefficient
    v = Lneg(2*k+1, chi12, f); v = mp.mpf(v.numerator)/v.denominator
    return v/((mp.mpf(f)/mp.pi)**(2*k+mp.mpf(3)/2)*mp.factorial(k)*mp.gamma(k+mp.mpf(3)/2)*(-1)**(k+1)/mp.pi)
K = 140
rem = [Lpos(k) for k in range(K)]
stokes = {}
for m in range(1, 14):
    c = rem[-1]*mp.mpf(m)**(2*(K-1)+2)
    stokes[m] = int(mp.nint(c)) if abs(c - mp.nint(c)) < 1e-6 else mp.nstr(c, 6)
    rem = [rem[k] - stokes[m]*mp.mpf(m)**(-(2*k+2)) if isinstance(stokes[m], int) else rem[k] for k in range(K)]
print('Stokes constants peeled from coefficients, m=1..13:', stokes, ' chi12:', [chi12(m) for m in range(1, 14)])
# (c) lowest zero of L(s,chi12) and the effect of moving it off the line
g1 = mp.findroot(lambda t: L(mp.mpf(1)/2 + 1j*t), 3.8)
print('lowest zero of L(s,chi12): 1/2 + i*', mp.nstr(mp.re(g1), 15))
rho = mp.mpf(1)/2 + 1j*mp.re(g1)
for d in [mp.mpf('0.01'), mp.mpf('0.1'), mp.mpf('0.3')]:
    rel = []
    for k in [0, 1, 2, 5, 10, 20]:
        s = 2*k + 2
        R = (1 - d**2/(s - rho)**2)*(1 - d**2/(s - mp.conj(rho))**2)
        rel.append(mp.nstr(mp.re(R) - 1, 4))
    print(f'delta={d}: relative change of T_k (k=0,1,2,5,10,20):', rel)
print('T_0 = 1 exactly (F(1)=1); moved zero with delta=0.1 would give T_0 =', mp.nstr(mp.re((1 - mp.mpf('0.01')/(2 - rho)**2)*(1 - mp.mpf('0.01')/(2 - mp.conj(rho))**2)), 15), '(not an integer)')
