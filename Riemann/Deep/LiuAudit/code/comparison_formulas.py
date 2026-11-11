#!/usr/bin/env python3
"""Exact checks of the explicit formulas of Liu 2610.26835v2 §13.4 (count crossing, l.9226-9269) and §14.1
(parameters, full exponent, endpoints, l.10016-10111) in number_fields_v4.tex."""
import sympy as sp
from fractions import Fraction as Fr
res = []
def rep(name, ok):
    res.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name, flush=True)
u, r, m, t, c, d, dl, x, b, ell, kap, a, g, th, z0, R = sp.symbols('u r m t c d delta x b ell kappa a g theta z0 R', real=True)
# §13.4 crossing: delta[u+(1-u)r] = delta[cu/2+(2-cu)(t-r)]
rs = sp.solve(sp.Eq(u+(1-u)*r, c*u/2+(2-c*u)*(t-r)), r)[0]
rep('l.9229 crossing r_* = [(2-cu)t+(c/2-1)u]/[3-(1+c)u]', sp.simplify(rs-((2-c*u)*t+(c/2-1)*u)/(3-(1+c)*u)) == 0)
ok = True; mn = Fr(10)
for ci in range(0, 41):
    cc = Fr(2, 3)+(Fr(12, 13)-Fr(2, 3))*ci/40
    for ui in range(0, 21):
        uu = Fr(ui, 40)
        for ti in range(0, 21):
            tt = 1+Fr(ti, 40)
            rr = ((2-cc*uu)*tt+(cc/2-1)*uu)/(3-(1+cc)*uu)
            mm = tt-rr
            mn = min(mn, rr)
            if not ((6-cc)/(10-2*cc) <= rr <= 1 and Fr(1, 3) <= mm <= Fr(1, 2) and rr >= Fr(8, 13)): ok = False
rep(f'l.9232-9234 (6-c)/(10-2c) <= r_* <= 1, 1/3 <= m_* <= 1/2, r_* >= 8/13 on c in [2/3,12/13], u in [0,1/2], t in [1,3/2] (min r_* = {mn})', ok and mn == Fr(8, 13))
ok = all(((4-cc)/(20-4*cc) <= Fr(5, 26) and 1-((6-cc)/(10-2*cc)) == (4-cc)/(10-2*cc) and 2*((6-cc)/(10-2*cc))-1 == 1/(5-cc) and cc/12 <= Fr(1, 13))
         for cc in [Fr(2, 3)+(Fr(12, 13)-Fr(2, 3))*k/50 for k in range(51)])
rep('l.9238-9249 z_M <= (4-c)/(20-4c) <= 5/26 < 1/5, 2r-1 >= 1/(5-c), z_S <= c/12 <= 1/13', ok)
rep('l.9248 c = 2/(3 kappa) in [2/3, 12/13] for kappa in [13/18, 1]', sp.Rational(2, 3)/sp.Rational(13, 18) == sp.Rational(12, 13))
# §14.1
lx = (1-b-ell)/2; ly = (1+b-ell)/2; h = (1+b+3*ell)/2; Bl = sp.Rational(11, 12)-ell/4; C = lambda s: s-(4+b)/6
rep('l.10018 h = 1-l_x+ell = l_y+2 ell', sp.simplify(h-(1-lx+ell)) == 0 and sp.simplify(h-(ly+2*ell)) == 0)
rep('l.10019 C(B(ell)) = l_x/2 + b/12', sp.simplify(C(Bl)-(lx/2+b/12)) == 0)
ck = 1/(3*kap)
D = 3-(1+2*ck)*x; P = (2-2*ck*x)*(1-x); Jc = (sp.Rational(5, 6)-dl)*D+dl*P; Rk = 1-dl+(sp.Rational(5, 6)-dl)*dl*P/(2*Jc)
rep('l.10050 2D-3P = 2x(2+c_k-3c_k x)', sp.simplify(sp.expand(2*D-3*P-2*x*(2+ck-3*ck*x))) == 0)
# with c=2c_k the count (v:cmp:10) matches (eq:v-refined-count)
cc_ = 2*ck
D2 = 3-(1+cc_)*x; P2 = (2-cc_*x)*(1-x)
rep('l.10044 (v:cmp:10) with c = 2c_kappa equals (eq:v-refined-count)', sp.simplify(D2-D) == 0 and sp.simplify(P2-P) == 0)
ok = True
for ki in range(0, 11):
    kk = Fr(13, 18)+(1-Fr(13, 18))*ki/10
    c_ = 1/(3*kk)
    for di in range(1, 51):
        dd = Fr(5, 6)*di/50
        for xi in range(0, 21):
            xx = Fr(xi, 40)
            DD = 3-(1+2*c_)*xx; PP = (2-2*c_*xx)*(1-xx); JJ = (Fr(5, 6)-dd)*DD+dd*PP
            RR = 1-dd+(Fr(5, 6)-dd)*dd*PP/(2*JJ)
            if not (DD > 0 and PP > 0 and Fr(1, 6) <= 1-dd <= RR <= 1-2*dd/3): ok = False
rep('l.10052 1/6 <= 1-delta <= R_kappa <= 1-2delta/3 on the grid (kappa, delta in (0,5/6], x in [0,1/2])', ok)
rep('l.10063 theta chain: (5/6)/(12*13/18) = 15/156 <= 1/10 and delta/(12 kappa) <= delta/2', sp.Rational(5, 6)/(12*sp.Rational(13, 18)) <= sp.Rational(1, 10) and 1/(12*sp.Rational(13, 18)) <= sp.Rational(1, 2))
# endpoint at d = h for delta >= 1/8: theta = c_k g/2 (= g/(6 kappa)), a = (1+delta)/2, g = x delta
Kc = (3*R+1)/4
E = a-Bl-ly*(sp.Rational(1, 2)+th)+h*(z0-sp.Rational(1, 6))-ell/2+ell*g+d*(Kc+th/2-z0)
Eh = E.subs(d, h).subs(th, ck*g/2).subs(a, (1+dl)/2).subs(g, x*dl)
C0 = -sp.Rational(1, 4)+b/6+5*ell/4
Eclaim = C0+dl/2+dl*x*(ell+ck*(h/4-ly/2))-3*h/4*(1-R)
rep('l.10091 endpoint E_{b,ell,kappa}(delta,x) at d=h (theta = c_k g/2 = g/(6kappa))', sp.simplify(sp.expand(Eh-Eclaim)) == 0)
rep('l.10091 theta = g/(6 kappa) equals c_kappa g/2', sp.simplify(ck*g/2-g/(6*kap)) == 0)
Es = E.subs(d, h).subs(th, g*ell/h).subs(a, (1+dl)/2).subs(g, x*dl)
Esc = C0+dl/2+ell*dl*x*(sp.Rational(1, 2)+2*ell/h)-3*h/4*(1-R)
rep('l.10100 small endpoint E_small at d=h (theta = g ell/h)', sp.simplify(sp.expand(Es-Esc)) == 0)
# (v:cmp:12) -> (eq:v-refined-full-exponent): U^{R+delta/2} -> U^{K+theta/2}, Y^{-a} -> Y^{-1/2-theta}
EB = a-sp.Symbol('B')+h*(z0-sp.Rational(1, 6))-a*ly-ell/2+g*ell+d*(R+dl/2-z0)
Efull = EB.subs(sp.Symbol('B'), Bl).subs(a*ly, 0)
E2 = a-Bl+h*(z0-sp.Rational(1, 6))-ly*(sp.Rational(1, 2)+th)-ell/2+g*ell+d*(Kc+th/2-z0)
rep('l.10076 full exponent = (v:cmp:12) with Y^{-a} -> Y^{-1/2-theta}, U^{R+delta/2} -> U^{K+theta/2}', sp.simplify(sp.expand((EB.subs(sp.Symbol('B'), Bl)+a*ly-ly*(sp.Rational(1, 2)+th)-d*(R+dl/2)+d*(Kc+th/2))-E2)) == 0)

# Thm 13.15 (l.9904-9918): B_inf = (18 xi+7)/(18 xi+9), xi the root of 1188 xi^3-2160 xi^2+171 xi+190 in [19429/50000, 38859/100000]
xs = sp.Symbol('xs')
pol = 1188*xs**3-2160*xs**2+171*xs+190
roots_in = [rt for rt in sp.Poly(pol, xs).intervals() if rt[0][0] <= sp.Rational(38859, 100000) and rt[0][1] >= sp.Rational(19429, 50000)]
sign_change = (pol.subs(xs, sp.Rational(19429, 50000))*pol.subs(xs, sp.Rational(38859, 100000))) < 0
xi_inf = sp.nsolve(pol, xs, 0.3886, prec=40)
Binf = (18*xi_inf+7)/(18*xi_inf+9)
rep(f'Thm 13.15 B_inf = {sp.N(Binf, 22)} (printed 0.874957019420098946; one sign change on the interval: {bool(sign_change)})', bool(sign_change) and abs(Binf-sp.Rational(874957019420098946, 10**18)) < 1e-18)


# Prop 14.3 (l.10261-10449): Hoelder block exponent, balance, maximal-length gate, omega
gg, vv, mm, LL, RR = sp.symbols('g v m L R', real=True)
thv = gg*vv
blk = (3*RR+1)/4 - thv*mm - gg*vv*(1-2*mm)/2
rep('l.10419 supply-limited: (3R+1)/4 - theta m - g v (1-2m)/2 = (3R+1)/4 - theta/2 (theta = g v)', sp.simplify(blk-((3*RR+1)/4-thv/2)) == 0)
k6 = 6*kap
m0 = sp.Rational(1, 2)-3*kap*LL
th2 = gg/(6*kap)
rep('l.10433 maximal-length: z_p(m0) = (1-2m0)/(6 kappa) = L', sp.simplify((1-2*m0)/(6*kap)-LL) == 0)
rep('l.10433 maximal-length: -theta m - theta(1-2m)/2 = -theta/2 for theta = g/(6kappa), z_p = (1-2m)/(6kappa)', sp.simplify(-th2*mm - gg*((1-2*mm)/(6*kap))/2 + th2/2) == 0)
gate = RR + (sp.Rational(1, 2)-th2)*m0 - ((3*RR+1)/4 - th2/2)
rep('l.10441 R + (1/2-theta) m0 - ((3R+1)/4 - theta/2) = (R - (6kappa-2g)L)/4', sp.simplify(sp.expand(gate - (RR-(6*kap-2*gg)*LL)/4)) == 0)
rep('l.10335 Hoelder: (3/4)R + (1-2g z_p)/4 = (3R+1)/4 - g z_p/2', sp.simplify(sp.Rational(3, 4)*RR+(1-2*gg*sp.Symbol('zp'))/4-((3*RR+1)/4-gg*sp.Symbol('zp')/2)) == 0)
epsv = sp.Symbol('epsv', positive=True)
om = 2*sp.Rational(5, 12)*epsv/100 + sp.Rational(5, 6)*epsv/100
rep('l.10329 omega = 2 g Delta + delta eta_M2 <= eps/60 (g <= 5/12, delta <= 5/6, Delta, eta_M2 <= eps/100)', sp.simplify(om - epsv/60) == 0)
rep('l.10404 sixth-power-free u: prod_{p|u} q_p >= q_u^(1/5) (exponents <= 5)', all(1 >= Fr(j, 5) for j in range(1, 6)))
rep('l.10313 2m + 6 kappa v(1-2m) <= 1 when v <= 1/(6 kappa)', sp.simplify((2*mm+6*kap*(1/(6*kap))*(1-2*mm)) - 1) == 0)


# Prop 14.6 (l.10672-10697)
val = sp.Rational(21, 100)*(6*sp.Rational(7497, 10000)-sp.Rational(1, 8))-(1-sp.Rational(2, 3)*sp.Rational(1, 8))
rep(f'l.10675 (21/100)(6 kappa - delta) - (1 - 2delta/3) at kappa = 7497/10000, delta = 1/8 equals {val} = 1279/750000, increasing in delta and kappa', val == sp.Rational(1279, 750000) and sp.Rational(2, 3)-sp.Rational(21, 100) > 0)
hh = sp.Symbol('hh', positive=True); C0s = sp.Symbol('C0s')
Esm_bound = C0s+dl/2+ell*dl*sp.Rational(1, 2)*(sp.Rational(1, 2)+2*ell/hh) - 3*hh/4*(2*dl/3)
Asm = sp.Rational(1, 2)+(ell/2)*(sp.Rational(1, 2)+2*ell/hh)-hh/2
rep('l.10688 E_small <= C0 + A_sm delta (x <= 1/2, 1-R >= 2delta/3)', sp.simplify(sp.expand(Esm_bound-(C0s+Asm*dl))) == 0)

print('SUMMARY_CMP', sum(res), 'pass of', len(res))

