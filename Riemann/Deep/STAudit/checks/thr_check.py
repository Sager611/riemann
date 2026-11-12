
import sympy as sp
e, y, v, nu, A, s = sp.symbols('eta y v nu A sigma', positive=True)
ok = []; bad = []
def chk(name, cond):
    (ok if cond else bad).append(name)
# m = 2 algebra
f = (44*e-1)/(16*e*(1+e))
crit = sp.solve(sp.diff(f, e), e)
crit = [c for c in crit if c.is_positive]
estar = sp.nsimplify((1+3*sp.sqrt(5))/44)
chk('critical point (1+3sqrt5)/44', len(crit)==1 and sp.simplify(crit[0]-estar)==0)
A2 = sp.nsimplify((23-3*sp.sqrt(5))/8)
chk('max f = (23-3sqrt5)/8', sp.simplify(f.subs(e, estar) - A2) == 0)
yy = sp.Rational(15,8)/(1+e)
F2 = 15*(5-4*s)/(16*(2-s))
chk('F2 = y(1/2+2eta)', sp.simplify((yy*(sp.Rational(1,2)+2*e)) - F2.subs(s, 1-e)) == 0)
chk('c2 = 15/8 - y(1/2-eta) = F2', sp.simplify(sp.Rational(15,8) - yy*(sp.Rational(1,2)-e) - F2.subs(s,1-e)) == 0)
chk('(F2-1)/eta = f', sp.simplify((F2.subs(s,1-e)-1)/e - f) == 0)
R = sp.Rational(3,2) - (yy/4)*(1-2*e)
sol = sp.solve(sp.Eq(R, F2.subs(s,1-e)), e)
chk('R(y/4) = F2 at eta = 1/14', sp.Rational(1,14) in sol)
chk('R(y/4) <= F2 for eta=0.2', (R - F2.subs(s,1-e)).subs(e, sp.Rational(1,5)) <= 0)
chk('R(y/4) > F2 for eta=0.05', (R - F2.subs(s,1-e)).subs(e, sp.Rational(1,20)) > 0)
chk('eta* in [1/14, 1/(2A2)]', sp.Rational(1,14) <= estar and estar <= 1/(2*A2))
chk('sigma0 = 1-1/(2A2) >= 3/4', 1 - 1/(2*A2) >= sp.Rational(3,4))
chk('A2 numeric 2.036474', abs(float(A2) - 2.0364745) < 1e-6)
# second branch of R and its left end
br2 = sp.Rational(3,4) + v + nu*(2*e - sp.Rational(1,2))
chk('left end of branch 2 = 2v+3eta-4v eta', sp.simplify(br2.subs(nu, sp.Rational(3,2)-2*v) - (2*v+3*e-4*v*e)) == 0)
chk('2v+3eta-4v eta <= 3/2-v(1-2eta) iff v<=1/2 (eta<1/2)', sp.simplify((sp.Rational(3,2)-v*(1-2*e)) - (2*v+3*e-4*v*e) - 3*(sp.Rational(1,2)-v)*(1-2*e)) == 0)
# lift branch at nu=y with v=y/4
chk('lift at nu=y: y/2+2y eta', sp.simplify(sp.Max(y/2, y/4+y/4) + 2*y*e - (y/2+2*y*e)) == 0)
# 13/14 margin
ee = sp.Rational(1,14)
chk('12e(13/14) = 3/4', sp.simplify(12*ee/(1+2*ee)) == sp.Rational(3,4))
# m = 1: class I and II give 2 - sigma; bar 1 + A(1-sigma)
chk('m=1: X^{3/2}V^{1-2s}, V=X^{1/2}: 2-s', sp.simplify(sp.Rational(3,2) + (1-2*s)/2 - (2-s)) == 0)
# m = 3
g = sp.Rational(3,4) + 6*e/(1+2*e)
chk('m=3 class II with y=3/(1+2eta) is 3/4+6e', sp.simplify(sp.Rational(9,4) - 3/(1+2*e)*(sp.Rational(1,2)-e) - g) == 0)
chk('m=3: 3/4+6e - A eta at eta=1/(2A) equals 1/4+3/(A+1)', sp.simplify(g.subs(e, 1/(2*A)) - A/(2*A) - (sp.Rational(1,4)+3/(A+1))) == 0)
chk('1/4+3/(A+1) < 1 iff A>3', sp.solve(sp.Eq(sp.Rational(1,4)+3/(A+1), 1), A) == [3])
# m = 4: y > 2A/(A-1) > 11/4 for A < 11/3; 3/4 + 11/(4A) < 3/2 iff A > 11/3
chk('2A/(A-1) at A=11/3 is 11/4', sp.simplify((2*A/(A-1)).subs(A, sp.Rational(11,3)) - sp.Rational(11,4)) == 0)
chk('3/4+11/(4A)=3/2 at A=11/3', sp.simplify((sp.Rational(3,4)+sp.Rational(11,4)/A).subs(A, sp.Rational(11,3)) - sp.Rational(3,2)) == 0)
# discrepancy: A > (2n+3)/3, n = ceil((3A-3)/2)-1 -> 5/(n+1) <= 10/(3A-3)
import math
okd = all(5/(math.ceil((3*a-3)/2)) <= 10/(3*a-3) + 1e-12 for a in [3.7+0.01*k for k in range(500)])
chk('discrepancy 5/(n+1) <= 10/(3A-3)', okd)
# A(J) bound: (M+3)/3 with M <= 75.4/L + 2  ->  <= 25.2/L + 1.7
L = sp.symbols('L', positive=True)
chk('(75.4/L+2+3)/3 <= 25.2/L + 1.7', sp.simplify(sp.Rational(252,10)/L + sp.Rational(17,10) - (sp.Rational(754,10)/L + 5)/3) .subs(L, 1) >= 0)
print('OK', len(ok)); print('FAIL', bad)
for n in ok: print('  ok:', n)
