#!/usr/bin/env python3
"""Exact checks of the explicit exponent identities/inequalities of Liu 2610.26835v2 §§12.5-12.8 and Thm 12.16
(number_fields_v4.tex). Each check prints PASS/FAIL with its TeX line."""
import sympy as sp
from fractions import Fraction as Fr
from itertools import product
res = []
def rep(name, ok):
    res.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name, flush=True)

M, A, z, q, xi, m, c, d, R, E, p, s0, K, K0, w, wo, Bc, Bd, g, ell, g2, t2, V, b2, p2, a0, J, d1, sig = sp.symbols(
    'M A z q xi m c d R E p s0 K K0 w wo Bc Bd g ell g2 t2 V b2 p2 a0 J delta1 sigma', real=True)
# §12.5 l.6283-6300: A_comp + q_k z <= (5/6+1/(3q_k))M + xi <= 14M/15 + xi when A > 5M/6, z <= (M-A)/q_k, q_k >= 10/3
qk = sp.symbols('q_k', positive=True)
expr = sp.Rational(3, 2)*M - A + (qk+2)*(M-A)/qk          # A_comp + q z at z = (M-A)/q (coefficient of z positive)
val = sp.simplify(expr.subs(A, sp.Rational(5, 6)*M))       # decreasing in A, so max at A = 5M/6
rep('l.6296 sup = (5/6+1/(3q))M: ' + str(sp.simplify(val/M)), sp.simplify(val - (sp.Rational(5, 6)+1/(3*qk))*M) == 0)
rep('l.6298 (5/6+1/(3q)) <= 14/15 for q >= 10/3', sp.Rational(5, 6)+sp.Rational(1, 3)/sp.Rational(10, 3) == sp.Rational(14, 15))
rep('l.6284 A_comp <= 23M/30: 3M/2-A+2z < 2M/3+M/10 at A>5M/6, z<M/20', sp.Rational(3, 2)-sp.Rational(5, 6)+sp.Rational(2, 20) == sp.Rational(23, 30))
rep('l.6282 6kappa-1 >= 10/3 <=> kappa >= 13/18', 6*sp.Rational(13, 18)-1 == sp.Rational(10, 3))
# §12.7.1 l.6736-6742 outer exponent identity
qt = q+R+E
lhs = (m-A-R/2-E)+(p+s0)+sp.Rational(1, 2)*(2*A-c-d+2*qt+Bc+Bd-2*s0)
rhs = (m+q)+sp.Rational(1, 2)*(Bc+Bd-(c+d-2*p-R))
rep('l.6736 complete outer exponent = M + (B_c+B_d-(c+d-2p-R))/2', sp.expand(lhs-rhs) == 0)
# l.6744-6750 brace nonpositive, per prime, i >= j >= 1, r = 1_{6 not| i-j}
ok = True
for i in range(1, 121):
    for j in range(1, i+1):
        r = 1 if (i-j) % 6 else 0
        pos = max(0, 3*i-5*j-r)
        neg_other = max(0, 3*j-5*i-r)
        if Fr(pos+neg_other, 6) > i+j-2-r: ok = False; print('  brace violated', i, j)
rep('l.6744 (3i-5j-r)_+/6 + (3j-5i-r)_+/6 <= i+j-2-r for 1<=j<=i<=120', ok)
# §12.7.2 diagonal excess l.6909-6915
c_expr = sp.solve(sp.Eq(K0, 2*A-c-d+R+E-m), c)[0]
Lam = a0+qt+wo+w+Bc-s0
diag = K+g-ell-s0
exc = sp.expand((diag-Lam).subs(a0, A-c-w).subs(c, c_expr))
target = sp.expand(A-(m+q)-d-(K0-K)-wo-Bc+g-ell)
rep('l.6909 diagonal excess = A-M-d-(K0-K)-w_o-B_c+g-ell', sp.simplify(exc-target) == 0)
# bound: -D0 - w_o + J_+ <= delta1 in both signs of J, with J = d-c+(K0-K)-2w+w_o, K <= K0+delta1, c,d,w>=0, w_o<=w
ok = True
import random
random.seed(1)
for _ in range(20000):
    cc, dd, ww = [Fr(random.randint(0, 300), 100) for _ in range(3)]
    wwo = Fr(random.randint(0, int(ww*100)), 100)
    dl = Fr(random.randint(0, 50), 100)
    KK0 = Fr(random.randint(0, 400), 100); KK = KK0+dl-Fr(random.randint(0, 300), 100)
    JJ = dd-cc+(KK0-KK)-2*ww+wwo
    if -(dd+KK0-KK)-wwo+max(JJ, 0) > dl: ok = False; break
rep('l.6913 -D0-w_o+J_+ <= delta1 (20000 random rational instances)', ok)
rep('l.6915 largest increment 2sigma - sigma/3 = 5sigma/3', 2*sig - sig/3 - sp.Rational(5, 3)*sig == 0)
# §12.7.3 child width l.6991-6995 and allowance l.7022-7030
Mp = (2*a0-K-g-g2-V)+(qt+wo+t2+V)
Jdef = d-c+(K0-K)-2*w+wo
e1 = sp.expand((Mp-((m+q)+Jdef-g-g2+t2)).subs(a0, A-c-w).subs(K0, 2*A-c-d+R+E-m))
rep('l.6993 M\' = m\'+q\' = M+J-g-g2+t2', sp.simplify(e1) == 0)
outside = (K+g-ell-a0)+(p2-s0)+(g2-t2)+(a0-b2)-a0
Dp = b2-p2+w+Bc+ell
rep('l.7022 allowance - outside exponent = M\' + Delta\'', sp.simplify(sp.expand(Lam-outside-(Mp+Dp))) == 0)
# M' <= M - sigma when g >= J_+ + sigma and g2 >= t2: J - g - g2 + t2 <= -sigma
ok = all((JJ - (max(JJ, 0)+s) - G2 + T2) <= -s for JJ in [Fr(k, 7) for k in range(-20, 21)] for s in [Fr(1, 10), Fr(1, 3)] for G2 in [Fr(0), Fr(1, 2), Fr(2)] for T2 in [Fr(0), Fr(1, 2), Fr(2)] if T2 <= G2)
rep('l.6995 M\' <= M - sigma (g = J_+ + sigma or + 2 sigma, t2 <= g2)', ok)
# local correlation bounds table l.6968-6975 from Lemma 12.10's closed form
Pp = sp.symbols('P', positive=True)
ok = True
for i in range(1, 25):
    for j in range(1, 25):
        cmin = min(i, j)
        if i == j:
            vals = {'unit': [Pp**(i-1)*(Pp-2) if i % 6 == 0 else Pp**(i-1)*(-1)], 'nonunit': [Pp**(i-1)*(Pp-1)]}
            bound = {'unit': Pp**(i-1) if i % 6 else Pp**i, 'nonunit': Pp**i}
            for key in vals:
                for v in vals[key]:
                    for PP in (13, 37, 1000):
                        if abs(v.subs(Pp, PP)) > bound[key].subs(Pp, PP): ok = False; print('  table', i, j, key)
        else:
            j0 = cmin
            v = Pp**(j0-1)*(Pp-1)*(1 if j0 % 6 == 0 else 0)
            for PP in (13, 37, 1000):
                if abs(v.subs(Pp, PP)) > PP**j0: ok = False
rep('l.6968 local |F| bounds P^{i-1}, P^i, P^i, P^{j0} implied by Lemma 12.10 (i,j<=24)', ok)
# §12.8 l.7106-7114 and l.7121-7126
rep('l.7110 amplifier error: w >= sigma/4 => w + sigma <= 5w', all(Fr(ww)+s <= 5*Fr(ww) for s in [Fr(1, 10)] for ww in [s/4, s/3, 7*s/3]))
rep('l.7104 main term g - J_+ = 2 sigma = 6 ell (ell = sigma/3)', 2*sig - 6*(sig/3) == 0)
kap, dz, eta, th = sp.symbols('kappa d_z eta theta_N', positive=True)
# removed length d_z <= defect/(6 kappa) + eta, defect <= 6(w+ell)+delta1+2theta; cost kappa d_z <= w+ell+delta1/6+theta/3+kappa*eta
cost = kap*((6*(w+ell)+d1+2*th)/(6*kap)+eta)
rep('l.7121 kappa d_z <= Delta\' + delta1/6 + theta_N/3 + eta (kappa<=1, Delta\'>=w+ell)', sp.simplify(sp.expand(cost-(w+ell+d1/6+th/3+kap*eta))) == 0)
# Thm 12.16 l.7647-7726 parameter arithmetic
eps, rho, delta, Cs, Dsh = sp.symbols('epsilon rho delta C_* D_sharp', positive=True)
# with xi, eps0, eta < eps/(512 C* D#): C*(xi+eps0) < eps/(256 D#), 2 D# C*(eta+xi+eps0) < 3 eps/256
rep('l.7713 C*(xi+eps0) < eps/(256 D#)', sp.simplify(Cs*2*eps/(512*Cs*Dsh) - eps/(256*Dsh)) == 0)
rep('l.7714 2D# C*(eta+xi+eps0) < 3eps/256', sp.simplify(2*Dsh*Cs*3*eps/(512*Cs*Dsh) - 3*eps/256) == 0)
rep('l.7711 eps/8 + eps/256 + 3eps/256 < eps', sp.Rational(1, 8)+sp.Rational(1, 256)+sp.Rational(3, 256) < 1)
rep('l.7656 tau + rho/6 = 7rho/6 + delta + 14sigma/3 + C*(..): 3+5/3 = 14/3', 3+sp.Rational(5, 3) == sp.Rational(14, 3))
# ranks
ok = True
for N in range(1, 30):
    D = 4*N+3
    for t in (0, 1):
        for b in range(0, N+1):
            for cc in (0, 1):
                r = 2*(N+1)*t+2*b+cc
                if r > D: ok = False
    # t=1 -> t=0
    if 2*(2*(0)+2*N+1)+1 > 4*N+3 or 2*(2*(N+1)+0+0)+0 < 4*N+4: ok = False
rep('l.7609-7643 rank r <= D = 4N+3; t-change target <= 4N+3 < 4N+4 <= source', ok)
# width decrease: M' <= M - sigma/2 lowers b(M)=ceil(4M/sigma) by >= 2
import math
ok = all(math.ceil(4*(Mv - sv/2)/sv) <= math.ceil(4*Mv/sv)-2 for sv in (0.1, 0.07) for Mv in [k*0.013 for k in range(1, 200)])
rep('l.7641 M\' <= M - sigma/2 => b(M\') <= b(M) - 2 (and R drops: 8 - 3 > 0)', ok)
rep('l.7629 z=0 centered call: M/4 + 5M/12 + xi = 2M/3 + xi < 5M/6 iff xi < M/6 (xi < rho/60 < M/60)', sp.Rational(1, 4)+sp.Rational(5, 12) == sp.Rational(2, 3))
rep('l.7655 base: A >= z, A+(6k-1)z <= M => 6 kappa z <= A+(6kappa-1)z <= M', all(6*k*zz <= Aa+(6*k-1)*zz for k in (Fr(13,18),Fr(3,4),Fr(1)) for zz in (Fr(0),Fr(1,10),Fr(1,3)) for Aa in (zz, zz+Fr(1,7), zz+Fr(2))))

# §12.8.2 l.7203-7262: exceptional excess identity, F_1 and F_2 bounds; §12.8.3 l.7316-7325 sup
f_, V_, w2 = sp.symbols('f V_ w2', real=True)
mp_ = 2*a0-K-g-g2-V_
qp_ = qt+wo+t2+V_
Dp_ = b2-p2+w+Bc+ell
exc_ex = (a0-b2) + (mp_-f_)/6 - (mp_+qp_) - Dp_
F1 = c/6+5*(d+K0-K)/6+w/3+qt/6+wo+Bc-5*g/6+ell
F2 = 2*b2-5*g2/6-p2+t2+V_/6+f_/6
claim = A-sp.Rational(5, 6)*(m+q)-F1-F2
e2 = sp.expand((exc_ex-claim).subs(a0, A-c-w).subs(K0, 2*A-c-d+R+E-m))
rep('l.7206 exceptional excess = A-5M/6-F_1-F_2 (+2theta+delta2/6+2eps0)', sp.simplify(e2) == 0)
Jx = d-c+(K0-K)-2*w+wo
rep('l.7220 F_1(g=J) = c+2w+q~/6+w_o/6+B_c+ell', sp.simplify(sp.expand(F1.subs(g, Jx)-(c+2*w+qt/6+wo/6+Bc+ell))) == 0)
ok = True
for _ in range(40000):
    cc, dd = Fr(random.randint(0, 300), 100), Fr(random.randint(0, 300), 100)
    sg = Fr(random.randint(1, 30), 1000); ww = Fr(random.randint(0, int(7*sg*1000/3)), 1000)
    wwo = Fr(random.randint(0, int(ww*1000)), 1000) if ww > 0 else Fr(0)
    dl = Fr(random.randint(0, 50), 10000)
    KK0 = Fr(random.randint(0, 400), 100); KK = KK0+dl-Fr(random.randint(0, 300), 100)
    RR = Fr(random.randint(0, 200), 100); qq = Fr(random.randint(0, 200), 100); EE = Fr(random.randint(0, 100), 100)
    qtt = qq+RR+EE
    BBc = max(Fr(0), (3*cc-5*dd-RR)/6)
    JJ = dd-cc+(KK0-KK)-2*ww+wwo
    lam = random.choice([Fr(0), sg/3])
    gg = max(JJ, 0) + random.choice([Fr(0), sg, 2*sg])
    FF1 = cc/6+5*(dd+KK0-KK)/6+ww/3+qtt/6+wwo+BBc-5*gg/6+lam
    if FF1 < Fr(2, 3)*(cc+ww)-3*sg-5*dl/6:
        ok = False; print('  F1 counterexample', cc, dd, ww, wwo, dl, KK0, KK, RR, qq, EE, gg); break
rep('l.7236 F_1 >= 2(c+w)/3 - 3 sigma - 5 delta1/6 (40000 random admissible rational instances)', ok)
ok = True
for i in range(1, 200):
    rows = [(Fr(7*i, 6), i)]
    if i % 6: rows.append((Fr(7*i-5, 6)+(Fr(1, 3) if i == 1 else 0), i))
    if i % 6 == 0: rows.append((Fr(7*i, 6)-1, i))
    for j0 in range(6, i, 6):
        rows.append((i+Fr(j0, 6)-1, Fr(i+j0, 2)))
    for F2v, bb in rows:
        if F2v < Fr(2, 3)*bb: ok = False
rep('l.7240-7254 local F_2 >= 2 b_2/3 for every row of the table, i < 200 (equality only i=1 nonunit)', ok)
ok = True
for i in range(1, 60):
    if i % 6:
        ok &= (2*i-Fr(5*i, 6)-1+1 == Fr(7*i, 6))
        ok &= (2*i-Fr(5*i, 6)-1+Fr(1, 6)+(Fr(2, 6) if i == 1 else 0) == Fr(7*i-5, 6)+(Fr(1, 3) if i == 1 else 0))
    else:
        ok &= (2*i-Fr(5*i, 6)-1 == Fr(7*i, 6)-1)
    for j0 in range(6, i, 6):
        ok &= ((i+j0)-Fr(5*j0, 6)-1 == i+Fr(j0, 6)-1)
rep('l.7240 local F_2 values follow from F_2 = 2b2-5g2/6-p2+t2+V/6+f/6 (V=1 nonunit, f=2 at i=1)', ok)
vals = [(-Fr(2, 3)*vv - max(Fr(0), Fr(1)-vv)) for vv in [Fr(k, 100) for k in range(0, 400)]]
rep('l.7322 sup_v {-2v/3-(L-v)_+} = -2L/3 (L=1 grid)', max(vals) == Fr(-2, 3))

print('SUMMARY_EXP', sum(res), 'pass of', len(res))

