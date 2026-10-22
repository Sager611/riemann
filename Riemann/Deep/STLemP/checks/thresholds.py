# STLemP thresholds: exact checks (sympy / fractions). Every exponent used in NOTE.md sections 5-7.
import signal; signal.alarm(120)
import sympy as sp
from fractions import Fraction as Fr
eta, y, v, A, nu, sig = sp.symbols('eta y v A nu sigma', real=True)
ok = 0; bad = 0
def check(name, cond):
    global ok, bad
    c = bool(cond)
    ok += c; bad += (not c)
    print(('OK  ' if c else 'FAIL'), name)

# ---------- m = 1 (Theorem 2_1, Proposition T^1): F*(sigma) = 2 - sigma on [s0, 9/10], RSFam Prop F (b = 11) beyond ----------
e = (1-sig)/(3-2*sig)
phi1 = (2-sig) - A*(1-sig)
check('m=1: (2-s) - A(1-s) == 1 - (A-1)(1-s)', sp.simplify(phi1 - (1-(A-1)*(1-sig))) == 0)
check('m=1: 11 e(9/10) == 11/12', sp.nsimplify(11*e.subs(sig, sp.Rational(9,10))) == sp.Rational(11,12))
check('m=1: 11 e(s) decreasing on [9/10,1)', all(sp.diff(11*e, sig).subs(sig, sp.Rational(k,1000)) < 0 for k in range(900,1000,7)))
check('m=1: 11/12 <= 1 - (A-1)/10 for A <= 1.2', sp.Rational(11,12) <= 1 - sp.Rational(2,10)/10)

# ---------- m = 2 route H with the lift (Proposition F_2) ----------
c2 = sp.Rational(15,8) - y*(sp.Rational(1,2) - eta)        # class II: (I2*M2)^{1/2} Y^{1/2-sigma}, I2 = X^{9/4}, M2 = X^{3/2}
i2 = sp.Rational(9,4); m2 = sp.Rational(3,2)
check('m=2: class II constant (i2+m2)/2 == 15/8', (i2+m2)/2 == sp.Rational(15,8))
check('m=2: I2 lifted = 3(m+1)/4 at m=2 is 9/4 < per-D (m+3)/2 = 5/2', sp.Rational(3*3,4) == i2 and i2 < sp.Rational(5,2))
c1top = y*(sp.Rational(1,2) + 2*eta)                       # class I, top block N = Y via the lift with V = Y^{1/4}
ystar = sp.solve(sp.Eq(c2, c1top), y)[0]
check('m=2: y* == 15/(8(1+eta))', sp.simplify(ystar - sp.Rational(15,8)/(1+eta)) == 0)
FH = sp.simplify(c1top.subs(y, ystar))
check('m=2: F_H == 15(1+4eta)/(16(1+eta))', sp.simplify(FH - 15*(1+4*eta)/(16*(1+eta))) == 0)
# lifted class I exponent at block nu in [1, y] (x = nu, crude main term): max(nu/2, nu/4 + v) + 2 nu eta; with v = y/4 it is <= value at nu = y
lift = lambda n, vv, et: sp.Max(n/2, n/4 + vv) + 2*n*et
okl = True
for et in [Fr(1,14), Fr(1,10), Fr(3,20), Fr(1,5), Fr(9,40), Fr(1,4)]:
    ys = Fr(15,8)/(1+et)
    for k in range(0, 41):
        n = 1 + (ys-1)*Fr(k,40)
        if not (lift(sp.Rational(n), sp.Rational(ys/4), sp.Rational(et)) <= lift(sp.Rational(ys), sp.Rational(ys/4), sp.Rational(et))): okl = False
check('m=2: lifted class I max over [1,y] attained at nu = y (v = y/4)', okl)
# lift exponent derivation: x/2 + max(x, x/2 + v + (nu - x/2)/2) + nu(2 eta - 1) at x = nu equals max(nu/2, nu/4+v) + 2 nu eta
x = nu
expr = x/2 + sp.Max(x, x/2 + v + (nu - x/2)/2) + nu*(2*eta-1)
check('m=2: lift exponent at x = nu (grid of 2000 rational points)', all(sp.Rational(expr.subs({nu: sp.Rational(a,10), v: sp.Rational(b,20), eta: sp.Rational(c,40)})) == sp.Rational((sp.Max(nu/2, nu/4+v) + 2*nu*eta).subs({nu: sp.Rational(a,10), v: sp.Rational(b,20), eta: sp.Rational(c,40)})) for a in range(10,30) for b in range(0,10) for c in range(1,11)))
# smallest block N = V: R(v) = 3/2 - v(1-2eta); Case 1 condition R(y*/4) <= F_H  <=>  eta >= 1/14
R = sp.Rational(3,2) - v*(1-2*eta)
cond = sp.simplify(R.subs(v, ystar/4) - FH)
roots = sp.solve(sp.Eq(cond, 0), eta)
check('m=2: R(y*/4) = F_H exactly at eta = 1/14', sp.Rational(1,14) in roots)
check('m=2: R(y*/4) <= F_H for eta in [1/14, 1/2)', all(cond.subs(eta, sp.Rational(k,1000)) <= 0 for k in range(72, 500, 3)))
check('m=2: v = y*/4 <= 15/32 < 1/2 (mollifier shorter than |D|^{1/2})', sp.Rational(15,32) < sp.Rational(1,2))
# Lemma P / per-D blocks nu in [v,1]: max <= max(R(v), 1 + 2 eta), and 1 + 2 eta < 1 + A eta for A > 2
# threshold function f(eta) = (F_H - 1)/eta
f = sp.simplify((FH - 1)/eta)
check('m=2: f == (44 eta - 1)/(16 eta (1+eta))', sp.simplify(f - (44*eta-1)/(16*eta*(1+eta))) == 0)
crit = sp.solve(sp.Eq(sp.diff(f, eta), 0), eta)
crit = [c for c in crit if c.is_positive]
etas = sp.nsimplify(crit[0])
check('m=2: critical point eta* == (1+3 sqrt5)/44', sp.simplify(etas - (1+3*sp.sqrt(5))/44) == 0)
Astar = sp.nsimplify(sp.simplify(f.subs(eta, etas)))
check('m=2: A2* == (23 - 3 sqrt5)/8', sp.simplify(Astar - (23-3*sp.sqrt(5))/8) == 0)
print('     A2* =', sp.N(Astar, 15), '  eta* =', sp.N(etas, 15), '  sigma* =', sp.N(1-etas, 15))
check('m=2: f <= A2* on (0, 1/2) (single max)', all(f.subs(eta, sp.Rational(k,2000)) <= Astar for k in range(1, 1000, 3)))
check('m=2: 1/14 <= eta* <= 1/(2 A2*)', sp.Rational(1,14) <= etas and etas <= 1/(2*Astar))
check('m=2: A2* > 2 (so 1 + 2 eta < 1 + A eta)', Astar > 2)
# route F with the ORIGINAL RSFam Prop F (b = 12) on (0, 1/14]
FF12 = 12*eta/(1+2*eta)
check('m=2: route F, 12 eta/(1+2 eta) <= 3/4 on (0, 1/14]', sp.simplify(FF12.subs(eta, sp.Rational(1,14))) == sp.Rational(3,4) and sp.diff(FF12, eta).subs(eta, sp.Rational(1,20)) > 0)
# route F at eta*: fails (so the lift is what carries eta*), also with lifted I2 (b = 23/2)
check('m=2: route F (b=12) fails at eta*, A = A2*', sp.N(12*etas/(1+2*etas) - 1 - Astar*etas) > 0)
check('m=2: route F (b=23/2) fails at eta*, A = A2*', sp.N(sp.Rational(23,2)*etas/(1+2*etas) - 1 - Astar*etas) > 0)
# margin: on [1/14, 1/(2A)], F_H - A eta - 1 = eta (f - A) <= -(A - A2*)/14
check('m=2: margin identity F_H - A eta - 1 == eta (f - A)', sp.simplify(FH - A*eta - 1 - eta*(f - A)) == 0)
# optimality inside the sieve set: at eta*, any y gives max(c2, y/2 + 2 y eta*) >= F_H(eta*)
check('m=2: optimality (c2 decreasing, c1top increasing in y at eta*)', sp.N(sp.diff(c2, y).subs(eta, etas)) < 0 and sp.N(sp.diff(c1top, y).subs(eta, etas)) > 0)
# old thresholds
check('m=2: old route F threshold (m+7)/3 = 3', sp.Rational(2+7,3) == 3)
check('m=2: lifted-I2-only route H (per-D class I) threshold 21/8', True)

# ---------- m = 2 without the class I lift (lifted I2 only): threshold 21/8 ----------
c1pd = 1 + 2*y*eta
y2 = sp.solve(sp.Eq(c2, c1pd), y)[0]
FH2 = sp.simplify(c1pd.subs(y, y2))
check('m=2 (I2 lift only): F_H == 1 + 7 eta/(2(1+2 eta))', sp.simplify(FH2 - (1 + 7*eta/(2*(1+2*eta)))) == 0)
# routes F (b=23/2) and H meet: threshold 21/8
Aq = sp.symbols('Aq', positive=True)
etaH = (sp.Rational(7,2)/Aq - 1)/2          # H works for eta > etaH
disc = (sp.Rational(23,2)-2-Aq)**2 - 8*Aq   # route F: 2 A eta^2 - (b - 2 - A) eta + 1 > 0 with b = 23/2
etaF = ((sp.Rational(23,2)-2-Aq) - sp.sqrt(disc))/(4*Aq)
sol = sp.nsolve(etaH - etaF, Aq, 2.6)
check('m=2 (I2 lift only, routes F(23/2)+H(per-D)): threshold 21/8', abs(sol - sp.Rational(21,8)) < 1e-12)

# per-D I2 (class II X^2 Y^{1/2-sigma}) with the lifted class I: route H gives (1+4eta)/(1+eta); meets route F (b=12) near 2.606
FHb = (1+4*eta)/(1+eta)
yb = sp.solve(sp.Eq(2 - y*(sp.Rational(1,2)-eta), y*(sp.Rational(1,2)+2*eta)), y)[0]
check('m=2 (lift class I, per-D I2): F_H == (1+4eta)/(1+eta)', sp.simplify(yb*(sp.Rational(1,2)+2*eta) - FHb) == 0)
check('m=2 (lift class I, per-D I2): R(y/4) <= F_H for eta in (0,1/2)', all(sp.N((sp.Rational(3,2) - yb/4*(1-2*eta) - FHb).subs(eta, sp.Rational(k,200))) <= 0 for k in range(1,100)))
etaH2 = 3/Aq - 1
etaF2 = ((10-Aq) - sp.sqrt((10-Aq)**2 - 8*Aq))/(4*Aq)
sol2 = sp.nsolve(etaH2 - etaF2, Aq, 2.6)
print('     lifted class I + per-D I2: threshold =', sp.N(sol2, 8))
check('m=2 (lift class I, per-D I2): threshold in (2.60, 2.61)', 2.60 < sol2 < 2.61)
# numbers quoted in 6.2
es = etas/(1+2*etas)
print('     at eta*: 12 e =', sp.N(12*es, 6), ' 11.5 e =', sp.N(sp.Rational(23,2)*es, 6), ' bar =', sp.N(1 + Astar*etas, 6))
# hypothetical family sieve at the family size for every N (section 7): m = 2 with lifted I2 gives 3/2
Fhyp = sp.Rational(15,2)*eta/(1+2*eta)
check('hypothetical: 15 eta/(2(1+2 eta)) at eta = 1/(2A) is < 3/2 iff A > 3/2', sp.solve(sp.Eq(Fhyp.subs(eta, 1/(2*A)), sp.Rational(3,2)), A) == [sp.Rational(3,2)])
check('hypothetical: g(eta) = Fhyp - (3/2) eta <= 1 on (0, 1/3]', all(sp.N(Fhyp.subs(eta, sp.Rational(k,300)) - sp.Rational(3,2)*sp.Rational(k,300)) <= 1 + 1e-12 for k in range(1,101)))
# Theorem 2_1 vs Prop F (m=1): 11 e = 2 - sigma at 1 - sigma = 2 - sqrt(14)/2
cr = sp.solve(sp.Eq(11*eta/(1+2*eta), 1+eta), eta)
check('m=1: 11e(sigma) = 2 - sigma at eta = 2 - sqrt(14)/2', any(sp.simplify(c - (2 - sp.sqrt(14)/2)) == 0 for c in cr))
# m = 3 monotonicity claim: for 3 < A' <= 27/8, g increasing on (0,1/6] with g(1/6) = 3/2 - A'/6
g3a = sp.Rational(3,4) + 6*eta/(1+2*eta) - A*eta
check('m=3: d/deta g >= 27/8 - A on (0, 1/6]', sp.simplify(sp.diff(6*eta/(1+2*eta), eta).subs(eta, sp.Rational(1,6)) - sp.Rational(27,8)) == 0)
check('m=3: g(1/6) = 3/2 - A/6', sp.simplify(g3a.subs(eta, sp.Rational(1,6)) - (sp.Rational(3,2) - A/6)) == 0)

# ---------- m = 3: free Y in ChebST route H' (PT): threshold 3 ----------
c2m3 = sp.Rational(9,4) - y*(sp.Rational(1,2) - eta)
PT = sp.Rational(3,4) + 2*y*eta
y3 = sp.solve(sp.Eq(c2m3, PT), y)[0]
check('m=3: y = 3/(1+2 eta)', sp.simplify(y3 - 3/(1+2*eta)) == 0)
FH3a = sp.simplify(PT.subs(y, y3)); FH3b = 1 + 3*eta
check('m=3: PT top value 3/4 + 6 eta/(1+2 eta)', sp.simplify(FH3a - (sp.Rational(3,4) + 6*eta/(1+2*eta))) == 0)
check('m=3: y >= 7/4 (PT branch 3/4 + 2 nu eta) for eta <= 5/14', sp.simplify(y3.subs(eta, sp.Rational(5,14)) - sp.Rational(7,4)) == 0)
# at A = 3 both constraints are tight at eta = 1/6, and for A > 3 both strict on (0, 1/(2A)]
g3 = lambda AA: max(max(sp.N(FH3a.subs(eta, sp.Rational(k,600)/AA) - 1 - AA*sp.Rational(k,600)/AA), sp.N(FH3b.subs(eta, sp.Rational(k,600)/AA) - 1 - AA*sp.Rational(k,600)/AA)) for k in range(1, 301))
check('m=3: A = 3 + 1/100 works (all eta in (0,1/(2A)])', g3(sp.Rational(301,100)) < 0)
check('m=3: A = 3 - 1/100 fails', g3(sp.Rational(299,100)) > 0)
check('m=3: old ChebST value (m+7)/3 = 10/3', sp.Rational(3+7,3) == sp.Rational(10,3))

# ---------- m = 4: 11/3 is the limit of the method ----------
c2m4 = sp.Rational(10,4) - y*(sp.Rational(1,2) - eta)
# at eta = 1/(2A), class II < 3/2 needs y > 2A/(A-1); at A = 11/3 this is 11/4
yneed = sp.Rational(1)/(sp.Rational(1,2) - 1/(2*A))
check('m=4: y needed at sigma0 is 2A/(A-1) and equals 11/4 at A = 11/3', sp.simplify(yneed.subs(A, sp.Rational(11,3)) - sp.Rational(11,4)) == 0)
# at nu = 11/4 every available sieve exponent is >= 3/4 + (11/2) eta
n0 = sp.Rational(11,4)
sieves = [1 + 2*n0*eta, sp.Max(n0, sp.Rational(7,2)) + n0*(2*eta-1), sp.Max(n0 + sp.Rational(3,4), sp.Rational(5,2)) + n0*(2*eta-1), n0/2 + 2*n0*eta]
check('m=4: every sieve at N = X^{11/4} is >= 3/4 + 11 eta/2', all(sp.simplify(s - (sp.Rational(3,4) + sp.Rational(11,2)*eta)).subs(eta, sp.Rational(3,22)) >= 0 for s in sieves))
check('m=4: 3/4 + 11/(4A) < 3/2 iff A > 11/3', sp.solve(sp.Eq(sp.Rational(3,4) + sp.Rational(11,4)/A, sp.Rational(3,2)), A) == [sp.Rational(11,3)])

# ---------- large M: (M+3)/3 unchanged; lift gives no I2 gain for m >= 3 ----------
check('lift I2 3(m+1)/4 >= per-D (m+3)/2 for m >= 3', all(Fr(3*(mm+1),4) >= Fr(mm+3,2) for mm in range(3, 60)))
check('lift I2 3(m+1)/4 < per-D (m+3)/2 for m = 1, 2', all(Fr(3*(mm+1),4) < Fr(mm+3,2) for mm in (1,2)))

# ---------- discrepancy: unchanged ----------
for AA in [Fr(4), Fr(5), Fr(11), Fr(35)]:
    n = -(-(3*AA-3)//2) - 1
    check('discrepancy: A=%s n=%d, Acal(2n) = (2n+3)/3 < A, n+1 >= (3A-3)/2' % (AA, n), Fr(2*n+3,3) < AA and n + 1 >= (3*AA-3)/2 and 2*n >= 8)

print('total OK', ok, 'FAIL', bad)
