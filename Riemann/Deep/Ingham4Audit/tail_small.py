# Task (4)/(5): R_an(300) for theta = 0.856 recomputed from DERIVATION sec. 5.2's formula with independently computed inputs
# (A_max, B_max from kln_core's TeX formulas and the same sigma-cell rule; M*, M1 from the kernel's definition by quadrature),
# the tail-lemma monotonicity conditions, and the exact small-x arithmetic. mpmath 256 bits.
import math
from mpmath import mp, mpf, log, exp, sin, pi, quad, findroot, e as E
mp.prec = 256
exec(open(__file__.replace('tail_small.py', 'kln_core.py')).read())
rows = [r for r in kln_rows(mpf('3e12'), mpf('0.2561882')) if r['ok']]
B = mpf('0.8733819666132076'); sig0 = mpf('0.515'); ds = mpf('0.0025')
cells = []; s = sig0
while s < B:
    t = min(s + ds, B)
    best = min([r for r in rows if s > r['smin']], key=lambda r: 2*r['c0'] + r['c1']*(s + t))
    cells.append((s, t, best)); s = t
Amax = max(exp(c[2]['c0'] + c[2]['c1']*c[1]) for c in cells); Bmax = max(c[2]['q0'] - c[2]['q1']*c[0] for c in cells)
print('sigma-cells %d; Amax = %s, Bmax = %s (cert.py log: 9.58339, 6.10455)' % (len(cells), mp.nstr(Amax, 8), mp.nstr(Bmax, 8)))
s1 = mpf('0.2173'); a = mpf(1)/10; b = 1 - a; dm = mpf('1e-9')
def S(y):
    if y == 0: return mpf(1)
    return max(sin(y)/y, s1) if y <= pi else min(s1, 1/y)
def G(l): return (S(a*l/2) + 2*dm*min(a, 1/l))*(S(b*l/2) + 2*dm*min(b, 1/l)) if l > 0 else (1 + 2*dm*a)*(1 + 2*dm*b)
C2 = (2/a + 2*dm)*(2/b + 2*dm)
ys = findroot(lambda y: sin(y)/y - s1, mpf('2.6'))
BRK = sorted(set([2*ys/a, 2*ys/b, 2*pi/a, 2*pi/b, 2/(a*s1), 2/(b*s1), 1/a, 1/b]))
def M0(p):   # G(0) + p int_1^oo G l^{p-1} dl
    pts = [mpf(1)] + [x for x in BRK if x > 1]
    tot = sum(quad(lambda l: G(l)*l**(p - 1), [u, v]) for u, v in zip(pts, pts[1:]))
    return G(0) + p*tot + p*C2*pts[-1]**(p - 2)/(2 - p)
th = mpf(856)/1000; L = mpf(300); m = 2
om = 1 - th; y = exp(-om*L); cst = 1/(1 - y/2); Dlow = om*L; Dmax = om*L + log(cst)
pmax = mpf(math.ceil(float((mpf(8)/3*(1 - sig0) + 4/Dlow)*200) - 1e-9))/200; p1 = mpf(math.ceil(float((1 + 1/Dlow)*200) - 1e-9))/200
Mst = M0(pmax); M1 = M0(p1)
kap = 1 - 8*om/3; Xi = L - mpf(8)/3*Dmax - 2*log(Dmax)
cond = dict(lam0_lt1=y*mpf('3e12') < 1, main=(1 - B)*(kap - 2/L) > 4/L, Bterm=2/L < 1 - B, Sterm=th - sig0 > 2/L, PP=th > 3/L, kap=kap - 2/L > 0, Xi=Xi > 0, ydmax=y < dm)
J = (L + 1)/log(2)
Ran = 2*(L*Amax*Mst*Dmax**3*exp(-(1 - B)*Xi)/Xi + Bmax*Mst*Dmax**2*exp(-(1 - B)*L) + exp(-(1 - sig0)*L)*cst*exp(om*L)*Dmax*M1) \
      + m*(L + 1)*(E*exp(-L/2)/4 + J*E*exp(-2*L/3)/9 + log(J)*cst*exp(-th*L)) + 2*exp(-3*L)
print('p_max = %s, p1 = %s, M* = %s, M1 = %s, Xi = %s' % (pmax, p1, mp.nstr(Mst, 10), mp.nstr(M1, 10), mp.nstr(Xi, 10)))
print('conditions:', cond, 'all:', all(cond.values()))
print('R_an(300) = %s (cert.py: 0.01419639346)' % mp.nstr(Ran, 10))
# R_an at larger L (non-increasing, as the lemma says)
def Ran_at(LL):
    LL = mpf(LL); yy = exp(-om*LL); Dm = om*LL + log(cst); X = LL - mpf(8)/3*Dm - 2*log(Dm); JJ = (LL + 1)/log(2)
    return 2*(LL*Amax*Mst*Dm**3*exp(-(1 - B)*X)/X + Bmax*Mst*Dm**2*exp(-(1 - B)*LL) + exp(-(1 - sig0)*LL)*cst*exp(om*LL)*Dm*M1) \
        + m*(LL + 1)*(E*exp(-LL/2)/4 + JJ*E*exp(-2*LL/3)/9 + log(JJ)*cst*exp(-th*LL)) + 2*exp(-3*LL)
vals = [Ran_at(v) for v in (300, 310, 350, 500, 1000, 3000, 10**4)]
print('R_an at L = 300, 310, 350, 500, 1000, 3000, 1e4:', [mp.nstr(v, 6) for v in vals], 'non-increasing:', all(u >= v for u, v in zip(vals, vals[1:])))
# small x, exact integers
print('1476^1000 <= 9999991^856:', 1476**1000 <= 9999991**856, '; 2 e^42 < 4e18:', 2*exp(42) < mpf('4e18'), ' e^42 =', mp.nstr(exp(42), 6))
print('worst record ratio log 4/log 7 =', mp.nstr(log(4)/log(7), 12), '; 7 + 7^0.856 - 11 =', mp.nstr(7 + mpf(7)**th - 11, 6), '(x = 7: next prime 11 needs 7^theta >= 4)')
