"""C3 (N6 audit): x-family (GORZ J^{d,0}, coefficients C(d,j) gamma(j)) vs t-family (Polya/Spark: degree-m Jensen
polynomial of E(t) = F(-t^2), P_m(t) = sum_k C(m,k) E^(k)(0) t^k).  Exact (Fractions + flint exact real-root counts).
(1) single pair F(x) = (1+x/w)(1+x/conj w), w = 5 - i: thresholds of both families;
(2) Schur-Szego/MSS-multiplicative intertwining P_{2d}(X=-t^2) = J^{d,0} (.)_d M_d, M_d = sum C(d,j) 4^j (d-1/2)_j X^j:
    M_d real-rooted with negative roots (d <= 60), inverse multiplier not real-rooted (d = 2..60)."""
from fractions import Fraction as Fr
from math import comb
import flint
def nreal(coeffs):  # exact rational coefficients, ascending; returns (#real roots, degree)
    den = 1
    for c in coeffs: den = den*c.denominator//__import__('math').gcd(den, c.denominator)
    p = flint.fmpz_poly([int(c*den) for c in coeffs])
    rts = p.complex_roots()
    return sum(m for r, m in rts if r.imag == 0), p.degree()
w = complex(5, -1); re_inv = Fr(5, 26); abs2 = Fr(26)          # 1/w = (5+i)/26, |w|^2 = 26
g = [Fr(1), 2*re_inv, 2/abs2]                                   # F = 1 + 2Re(1/w) x + x^2/|w|^2 = sum g_j x^j/j!
xs = [d for d in range(2, 200) if nreal([comb(d, j)*g[j] for j in range(3)])[0] == 2]
print('x-family J^{d,0} real-rooted for d in [%d, %d] (gap-free: %s)' % (min(xs), max(xs), xs == list(range(2, max(xs)+1))))
# t-family: E(t) = F(-t^2) = 1 - 2Re(1/w) t^2 + t^4/|w|^2 ; E^(2)(0) = -2*2Re(1/w), E^(4)(0) = 24/|w|^2
E = {0: Fr(1), 2: -2*2*re_inv, 4: Fr(24)/abs2}
ts = []
for m in range(4, 400):
    c = [Fr(0)]*(m+1)
    for k, v in E.items(): c[k] = comb(m, k)*v
    r, deg = nreal(c)
    if r == deg: ts.append(m)
print('t-family P_m real-rooted for m in [%d, %d] (gap-free: %s)' % (min(ts), max(ts), ts == list(range(4, max(ts)+1))))
print('Spark criterion sin^2(arg w) =', Fr(1, 26), '<= (4m-6)/(m(m-1)) at m=103:', Fr(1,26) <= Fr(4*103-6, 103*102), ' m=104:', Fr(1,26) <= Fr(4*104-6, 104*103))
d = 27
c2d = [Fr(0)]*(2*d+1)
for k, v in E.items(): c2d[k] = comb(2*d, k)*v
print('d=27: J^{27,0} real roots', nreal([comb(d, j)*g[j] for j in range(3)]), ' P_54 real roots', nreal(c2d))
# (2) intertwining multiplier
def ff(a, j):
    r = Fr(1)
    for i in range(j): r *= (a - i)
    return r
bad_M, bad_inv = [], []
for d in range(1, 61):
    M = [comb(d, j)*4**j*ff(Fr(2*d-1, 2), j) for j in range(d+1)]
    r, deg = nreal(M)
    if r != deg: bad_M.append(d)
    Minv = [comb(d, j)/(4**j*ff(Fr(2*d-1, 2), j)) for j in range(d+1)]
    r2, deg2 = nreal(Minv)
    if r2 == deg2: bad_inv.append(d)
print('M_d not real-rooted for d in', bad_M, ' (coefficients > 0, so roots < 0)')
print('inverse multiplier real-rooted for d in', bad_inv, ' (2..60 checked)')
# identity check on a generic sequence: coefficient of X^j in P_{2d}: C(2d,2j) E^(2j)(0)(-1)^j with E^(2j)(0)=(-1)^j(2j)! g_j/j!
d = 7; gg = [Fr(k*k+1, 3**k) for k in range(d+1)]
lhs = [comb(2*d, 2*j)*Fr(__import__('math').factorial(2*j), __import__('math').factorial(j))*gg[j] for j in range(d+1)]
rhs = [comb(d, j)*gg[j]*4**j*ff(Fr(2*d-1, 2), j) for j in range(d+1)]
print('identity P_{2d} = J^{d,0} (.) M_d (d=7, generic sequence):', lhs == rhs)
