# b3pair range check 3 (second referee, 20 Nov 2026): R13 with B3N REPORT l.147 (§5.1 row 19) D > (1 + 2nu - theta + delta)/(2 eta), and R5's 0 < nu. Exact rationals.
from fractions import Fraction as F
k, l = F(89, 2742), F(1184, 1371); thc = F(4939, 10000)
nu1 = lambda th: (th*(1 + 2*l) - (k + l))/(2*(1 + l - k)); nu0 = nu1(thc)
print('nu0 =', nu0, '; R5 lower bound: nu0 > 0:', nu0 > 0)
for gap in (F(1, 10**4), F(1, 10**6)):
    th = thc + gap; eta = gap/2
    Dnew = (1 + 2*nu0 - th)/(2*eta); Dold = (1 + nu0 - th)/(2*eta)
    print(f'theta - theta_c = {gap}: D_min (B3N l.147, delta = 0) = {Dnew} = {float(Dnew):.6f}; old numerator (B3F l.110) = {Dold} = {float(Dold):.6f}; nu1(theta\') - nu0 = {nu1(th - eta) - nu0} > 0: {nu1(th - eta) - nu0 > 0}')
print('numerator 1 + 2nu0 - theta > 0 on (theta_c, 1):', 1 + 2*nu0 - thc > 0, 1 + 2*nu0 - 1 > 0, '; D_min finite for every fixed theta > theta_c')
