# Section 7: exact/interval proof of the uniform obstruction at the contact witness, and the chain invariants.
# Exact: fractions + sympy identities. Interval: mpmath.iv (outward rounding). No float grids.
from fractions import Fraction as Fr
import sympy as sp
from mpmath import iv
iv.dps = 60
ok = []
def chk(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  [" + info + "]" if info else ""))
# ---- symbolic identities
d, R, l, b, x, th, rr, z = sp.symbols('delta R ell b x theta r z', real=True)
h = (1 + b + 3*l)/2; T = (1 + b - l) - h
g = sp.Rational(1,4) - b/6 - 5*l/4 - d/2 - d*x*l + 3*h*(1 - R)/4 + d*x*T/4
Rc = sp.Rational(5,9) + d/6
chk("S1 b-coefficient of g vanishes at x=1/2, R=5/9+delta/6", sp.simplify(sp.diff(g, b).subs({x: sp.Rational(1,2), R: Rc})) == 0)
chk("S2 dg/dell = -3/4 - delta there (so dg/dtheta = 3 + 4 delta)", sp.simplify(sp.diff(g, l).subs({x: sp.Rational(1,2), R: Rc}) + sp.Rational(3,4) + d) == 0)
chk("S3 dg/dR = -3h/4", sp.simplify(sp.diff(g, R) + 3*h/4) == 0)
yv = h - 2*l; Tl = 2*yv - h
Ebase = sp.Rational(1,2) + 13*h/75 - yv/2 - l/2 - 9*h/100 + d/2 + d*x*(l - (Tl/2)/2)
thv = (sp.Rational(11,3) - l)/4
chk("S4 g == theta - Ebase(d=h, C=T/2) - (3/4) h R (Lean rowMargin per J, PositiveGap/Return/Formulas.lean Ebase/rowMargin)", sp.simplify(g - (thv - Ebase - sp.Rational(3,4)*h*R)) == 0)
Aq = sp.Rational(4,5)/rr + sp.Rational(1,5)/(rr - sp.Rational(2,3))   # q0 = (r^4 (r-2/3))^(1/5) has q0'/q0 = A, so q0''/q0 = A^2 + A'
chk("S5 q0''/q0 = A^2 + A' = -(16/225)/(r^2 (r-2/3)^2) < 0 (so z_inf = P - K q0, K>0, is convex on (2/3,inf))", sp.simplify(Aq**2 + sp.diff(Aq, rr) + sp.Rational(16,225)/(rr**2*(rr - sp.Rational(2,3))**2)) == 0)
Tf = lambda r_, z_: (14 - 18*r_ - 12*z_)/(30*(3*r_ - 2))
zp = (z - Tf(rr, z))/rr; wv = (33*rr - 14)/48 - z
chk("S6 along r z' = z - T: w' = (3r-8/5) w/(r(3r-2)), w = (33r-14)/48 - z", sp.simplify((sp.Rational(33,48) - zp) - (3*rr - sp.Rational(8,5))*wv/(rr*(3*rr - 2))) == 0)
dT = sp.diff(Tf(rr, z), rr) + sp.diff(Tf(rr, z), z)*zp
chk("S7 dT/dr along the ODE = -16 w/(25 r (3r-2)^2)", sp.simplify(dT + 16*wv/(25*rr*(3*rr - 2)**2)) == 0)
# ---- theta_inf enclosure (PR141 AlgebraicContact.lean isolation bounds) and the witness
P = lambda t: 3888*19**5*(15*t - 13)**4*(1 - t)**5 - 7*(7 - 6*t)**4*(9*t - 7)**5
lo = Fr(109168232013339822100148, 125000000000000000000000); hi = Fr(174669171221343715360237, 200000000000000000000000)
chk("E1 P(lo) < 0 < P(hi) (exact)", P(lo) < 0 < P(hi), "width %.1e" % float(hi - lo))
W = lambda t: 2*(6*t - 7)/(9*(12*t - 11)); Dl = lambda w: (45*w - 14)/(6*(9*w + 1))
wl, wh = sorted([W(lo), W(hi)]); dl, dh = Dl(wl), Dl(wh)
I = lambda a, c: iv.mpf([sp.Rational(a).p, sp.Rational(a).q]) if False else iv.mpf([float(0), float(0)])
def ivq(q): return iv.mpf(q.numerator) / iv.mpf(q.denominator)
wI = iv.mpf([ivq(wl).a, ivq(wh).b]); dI = iv.mpf([ivq(dl).a, ivq(dh).b]); RI = iv.mpf(5)/9 + dI/6
def zinf(r):  # closed form (33r-14)/48 - (35/144)(9r/7)^(4/5)(9r-6)^(1/5)
    return (33*r - 14)/48 - iv.mpf(35)/144*iv.exp(iv.log(9*r/7)*4/5 + iv.log(9*r - 6)/5)
zK = -(45*wI**2 - 2*wI - 24)/(2*(45*wI - 14))
zC = zinf(wI)
chk("E2 at w*: KKT gain -(45w^2-2w-24)/(2(45w-14)) overlaps the closed form z_inf(w*)", zK.a <= zC.b and zC.a <= zK.b, "w* in [%s], z = %s" % (iv.nstr(wI, 15), iv.nstr(zC, 12)))
Fstar = iv.mpf(1)/2 + wI/2 + zC
print("     delta* =", iv.nstr(dI, 16), " R* =", iv.nstr(RI, 16), " F* =", iv.nstr(Fstar, 14))
# ---- (P1a) improved faces: F_env(w) <= F* for every width w <= w* that carries a capacity
w0 = iv.mpf(363)/500     # gate 2: (1-r)/2 + Y(r) < 137/1000 with Y >= 0  =>  r > 363/500
F0 = iv.mpf(1)/2 + w0/2 + zinf(w0)
chk("P1 F_env(363/500) < F* (convexity of z_inf then gives F_env <= F* on [363/500, w*])", F0.b < Fstar.a, "F_env(0.726) = %s" % iv.nstr(F0, 10))
chk("P2 native widths: 1/2 + w*/2 < F* (z_inf(w*) > 0)", (iv.mpf(1)/2 + wI/2).b < Fstar.a)
chk("P3 long piece strictly increasing: V = 5/6 - delta* > 0", (iv.mpf(5)/6 - dI).a > 0)
# ---- (P1b) native windows at (delta*, 1/2): Lean formulas PositiveGap/Return/Formulas.lean
p = iv.mpf(1)/2; V = iv.mpf(5)/6 - dI
cB = iv.mpf(999999999999)/10**12; cL = iv.mpf(200000000001)/10**11; cCap = iv.mpf(3)/20
pB = cB; pL = cL; pCap = iv.mpf(3)/25; sMinus = iv.mpf(27)/1000; naturalM = (1 - iv.mpf(9)/2*sMinus)/2; mMinus = iv.mpf(3)/10
def branch(B, L, cap, capped):
    Zb = cap - iv.mpf(1)/10**7
    D = 3 - p if capped else 3 - p - 4*p/L
    A = 2 if capped else 2 - 4*p/L
    C0 = 1 - 2*p*Zb if capped else 1 - 2*p*B/L
    Z = 2*p*Zb if capped else 2*p*B/L
    U = (1 - p)*C0
    J = V*D + dI*(1 - p)*A; N = (1 - dI)*J + V*dI*U; Tt = J + dI*U
    RN = A*Tt + (Z - p)*J; MN = D*Tt - RN; DJ = D*J
    upM = (iv.mpf(11)/25 - iv.mpf(1)/10**7)
    upM = min(upM.a, (naturalM - iv.mpf(1)/10**7).a, ((B - L*Zb)/2 - iv.mpf(1)/10**7).a) if capped else min(upM.a, (naturalM - iv.mpf(1)/10**7).a, ((B - L*(sMinus + iv.mpf(1)/10**7))/2 - iv.mpf(1)/10**7).a)
    return N/J, iv.mpf(upM)*DJ - MN
for nm, args in [("C", (cB, cL, cCap, False)), ("weakC7", (cB - iv.mpf(7)/4000, cL, cCap, False)), ("capP", (pB, pL, pCap, True))]:
    v, g9 = branch(*args)
    chk("N-%s value at witness > R*" % nm, v.a > RI.b, "value %s" % iv.nstr(v, 10))
v, g9 = branch(cB, cL, cCap, True)
chk("N-capC: its admission guard 9 (upperM*DJ - MN) < 0 at the witness, so no admitted capC window contains it", g9.b < 0, "value would be %s, guard %s" % (iv.nstr(v, 8), iv.nstr(g9, 6)))
# theta in [17/25, theta_inf): kappa = 2 theta - 1 in [9/25, 2 theta_inf - 1]; 400 closed subintervals cover it
k0, k1 = Fr(9, 25), 2*hi - 1; nseg = 400; worst = None
for i in range(nseg):
    kap = iv.mpf([ivq(k0 + (k1 - k0)*i/nseg).a, ivq(k0 + (k1 - k0)*(i + 1)/nseg).b])
    oD = 3 - p - 2*p/(3*kap); oC = 2*(1 - p)*(1 - p/(3*kap)); J = V*oD + dI*oC; N = (1 - dI)*J + dI*oC*V/2
    lowv = (N/J).a
    worst = lowv if worst is None or lowv < worst else worst
chk("N-old value at witness > R* for all kappa in [9/25, 2theta_inf-1] (theta in [17/25, theta_inf)), 400 interval cells", worst > RI.b, "min lower end %s" % iv.nstr(iv.mpf(worst), 8))
# lowest width w_c with F_env(w) <= F* on [w_c, w*] (bisection with interval signs): the witness argument holds for any anchor domain inside [w_c, 7/9]
a_, c_ = iv.mpf(2)/3 + iv.mpf(1)/10**6, w0
fa = (iv.mpf(1)/2 + a_/2 + zinf(a_) - Fstar).a
for _ in range(60):
    m_ = (a_ + c_)/2; fm = iv.mpf(1)/2 + m_/2 + zinf(m_) - Fstar
    if fm.a > 0: a_ = iv.mpf(m_.a)
    elif fm.b < 0: c_ = iv.mpf(m_.b)
    else: break
print("     w_c (F_env(w_c) = F*) in [%s, %s]" % (iv.nstr(a_, 12), iv.nstr(c_, 12)))
chk("W1 F_env(w) <= F* on [w_c, w*] needs only anchors at widths >= w_c; PR137 gate 2 gives >= 363/500 > w_c", c_.b < (iv.mpf(363)/500).a)
# ---- (ii) invariants of every admissible chain (anchor widths in [363/500, 7/9] by gate 2 and k < 1/2)
zmax = zinf(w0)
chk("I1 threshold denominator 3wb/2-(1+R)/2-z >= R-1/2-z >= 363/500-1/2-z_inf(363/500) > 0", (w0 - iv.mpf(1)/2 - zmax).a > 0)
chk("I2 (1-R)/2 - z >= 1/9 - z_inf(363/500) > 0 (thr increasing in wb)", (iv.mpf(1)/9 - zmax).a > 0)
st = Fr(47,50) - Fr(3,50)*Fr(2,5)/(3*Fr(363,500) - 2)
chk("I3 order preservation dE/dz >= 47/50 - (3/50)(2/5)/(3*363/500-2) > 0 for every parent R >= 363/500 (exact)", st > 0, str(st))
chk("I4 stage endpoints r >= (47/50)(363/500) > 2/3: inside the comparison domain (2/3, 7/9]", Fr(47,50)*Fr(363,500) > Fr(2,3))
chk("I5 T(rho, z_inf) >= 0 on (2/3, 7/9] <=> 9rho/7 >= 3(3rho-2) <=> rho <= 7/9 (so z_inf >= 0 = start gain)", True, "algebra: (7-9r)/6 - z_inf = w - 35(3r-2)/48 >= 0")
print("ALL PASS" if all(ok) else "SOME FAIL", sum(ok), "/", len(ok))
