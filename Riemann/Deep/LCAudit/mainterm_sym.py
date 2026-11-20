#!/usr/bin/env python3
"""lcaudit: Wu's main term at the shift, term by term, symbolically (sympy).
Leading order of Wu Prop 3.1 (1802.09704 l.400-401, exponent (h,k)^{1+a+b} corrected) for a(n) = mu(n) P(log(y/n)/log y),
alpha = A/L, beta = B/L, log y = theta L, chi_0 mod q.  Run: uv run --no-project --with sympy python3 mainterm_sym.py
"""
import sympy as sp
x, t, s, u, th, R, L, A, B, X = sp.symbols('x t s u theta R L A B X', real=True)
a2, a3, a4 = sp.symbols('a2 a3 a4', real=True)
phq = sp.symbols('phi_q_over_q', positive=True)       # phi(q)/q
P = x + a2*(x**2 - x) + a3*(x**3 - x) + a4*(x**4 - x)  # general quartic with P(0)=0, P(1)=1
coef = sp.Poly(P, x).all_coeffs()[::-1]                # p_0..p_4
dP = sp.diff(P, x)
out = []
def say(*m):
    s_ = ' '.join(str(z) for z in m); print(s_); out.append(s_)

# Step 1: residue  (1/2 pi i) int e^{sX} (beta+s)/s^{i+1} ds = beta X^i/i! + X^{i-1}/(i-1)!   (checked by series)
for i in range(0, 5):
    bet = sp.Symbol('bet')
    expr = sp.exp(s*X)*(bet + s)/s**(i+1)
    res = sp.residue(expr, s, 0)
    want = bet*X**i/sp.factorial(i) + (X**(i-1)/sp.factorial(i-1) if i >= 1 else 0)
    assert sp.simplify(res - want) == 0
say("step 1 residue formula: OK for i = 0..4")

# Step 2: Dirichlet series of the h,k sum: prod_{p not | q} (1 - p^{-1-beta-s} - p^{-1-alpha-u} + p^{-1-s-u})
#   = [zeta(1+s+u)/(zeta(1+beta+s) zeta(1+alpha+u))] * A(...) * prod_{p|q} (1-p^{-1-s-u})^{-1}... ; leading order near 0:
#   (q/phi(q)) (beta+s)(alpha+u)/(s+u)   [A(0)=1 per prime p not | q checked symbolically below; p | q factor -> (1-1/p)^{-1}]
pp = sp.Symbol('p', positive=True)
xx, yy, zz = 1/pp, 1/pp, 1/pp
assert sp.simplify((1 - xx - yy + zz)*(1 - zz)/((1 - xx)*(1 - yy)) - 1) == 0
say("step 2 Euler factor A_p(0) = 1 for p not | q: OK; p | q factors give (q/phi(q)) at leading order")

# Step 3: y^{s+u}/(s+u) -> log y * int_0^1 y^{(s+u)t} dt; separate s and u integrals with Step 1, X = t log y:
ly = th*L
def side(shift):   # sum_i p_i i!/(log y)^i [shift X^i/i! + X^{i-1}/(i-1)!],  X = t*ly
    tot = 0
    for i, p_i in enumerate(coef):
        Xv = t*ly
        term = shift*Xv**i/sp.factorial(i) + (Xv**(i-1)/sp.factorial(i-1) if i >= 1 else 0)
        tot += p_i*sp.factorial(i)/ly**i*term
    return sp.expand(tot)
def S_lead(al, be):   # leading order of sum_{h,k} chi0 mu mu P P (h,k)^{1+al+be} h^{-1-be} k^{-1-al}, without the q/phi factor
    return sp.integrate(sp.expand(ly*side(be)*side(al)), (t, 0, 1))
S_ab = S_lead(A/L, B/L)
target = sp.integrate(sp.expand((dP + B*th*P)*(dP + A*th*P)).subs(x, t), (t, 0, 1))/(th*L)
assert sp.simplify(S_ab - target) == 0
say("step 3 S(alpha,beta) ~ (q/phi) (1/(theta L)) int (P'+B theta P)(P'+A theta P): OK (general quartic P)")

# Step 4: Wu's two terms at leading order (chi_0): L(1+z,chi_0) ~ (phi/q)/z, (qw/2pi)^{-alpha-beta} -> e^{-(A+B)}
T1 = (phq*L/(A + B)) * (S_lead(A/L, B/L)/phq)
T2 = sp.exp(-(A + B)) * (phq*(-L)/(A + B)) * (S_lead(-B/L, -A/L)/phq)   # second term: (h,k)^{1-a-b} h^{-1+a} k^{-1+b} = S(-beta,-alpha)
c = sp.simplify(T1 + T2)
I = lambda e: sp.integrate(sp.expand(e), (x, 0, 1))
c_formula = (I((dP + A*th*P)*(dP + B*th*P)) - sp.exp(-(A + B))*I((dP - A*th*P)*(dP - B*th*P)))/((A + B)*th)
assert sp.simplify(c - c_formula) == 0
say("step 4 q/phi(q) cancels in both terms; c(A,B) = [int(P'+A th P)(P'+B th P) - e^{-(A+B)} int(P'-A th P)(P'-B th P)]/((A+B) th): OK")

c1 = 1 + (sp.exp(2*R) - 1)/(2*R*th)*I((dP + th*R*P)**2)           # Conrey's c_1 (BCY l.338-345), Q = 1
cRR = c_formula.subs({A: R, B: R}); cmRmR = c_formula.subs({A: -R, B: -R})
assert sp.simplify(cmRmR - c1) == 0
say("step 5 left shift a = b = -R: c = c_1 (Conrey): OK")
assert sp.simplify(cRR - sp.exp(-2*R)*c1) == 0
say("step 6 right shift a = b = +R: c = e^{-2R} c_1 =: c_+ : OK  (this is the hinge input)")
sym = sp.simplify(c_formula - sp.exp(-(A + B))*c_formula.subs({A: -B, B: -A}, simultaneous=True))
assert sym == 0
say("step 7 functional-equation symmetry c(A,B) = e^{-(A+B)} c(-B,-A): OK")
lim0 = sp.limit(cRR.subs({a2: 0, a3: 0, a4: 0}), R, 0)
assert sp.simplify(lim0 - (1 + 1/th)) == 0
say("step 8 P = x, R -> 0: c -> 1 + 1/theta (Wu / Lamzouri LemmaWu l.1004): OK")
# Step 9: inf over P of int (P' + bP)^2, b = theta R: EL minimiser sinh(bx)/sinh(b), value b(1+coth b); c_+ at it = c*(R)
b = sp.Symbol('b', positive=True)
Ps = sp.sinh(b*x)/sp.sinh(b)
val = sp.simplify(sp.integrate(sp.expand((sp.diff(Ps, x) + b*Ps)**2).rewrite(sp.exp), (x, 0, 1)))
assert sp.simplify((val - b*(1 + sp.coth(b))).rewrite(sp.exp)) == 0
cstar_sym = sp.exp(-2*R) + (1 - sp.exp(-2*R))/(2*R*th)*(th*R)*(1 + sp.coth(th*R))
alt = 1 + (1 - sp.exp(-2*R))/(sp.exp(2*th*R) - 1)
assert sp.simplify((cstar_sym - alt).rewrite(sp.exp)) == 0
say("step 9 inf_P c_+ = e^{-2R} + (1-e^{-2R})(1+coth(theta R))/2 = 1 + (1-e^{-2R})/(e^{2 theta R}-1): OK")
open('logs/mainterm_sym.log', 'w').write('\n'.join(out) + '\n')
