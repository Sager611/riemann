"""liuaudit2: exponent steps of Lemma 13.3 (number_fields_v4.tex l.8238-8502), mean squares (13.6) = v:cmp:2, (13.7) = v:cmp:3.
(13.6) final step (l.8344-8362):
  E(e) = O/2 + H - S0 - B0 + z_a - u - l_b - 2e/3 - (T_d - v - 3 l_b - e)_+/2,  u = min(v, z_a, (v+z_a)/3),
  T_d = 2H + 2A0 + 2z_a - N_* - N0 - 3B0, N_* = 1 + l', l' = l - d, M' = 1 - l - 2d, Delta_H = M' - O - H >= 0,
  constraints 2A0 <= O, N0 <= A0, all lengths >= 0, z_a <= l'.
  Claims: identity O/2+Delta_H+S0+B0+T_d/4 = (2M'+2z_a-1-l'+2Delta_H)/4 + (2A0-N0+B0+4S0)/4;
          E(e) <= E(0) - e/6;  E(0) <= M' + (d+5l-1)_+/4   (giving v:cmp:2).
(13.7) (l.8376-8500): zero frequency Q'/Y'; nonprincipal sum Q'/Y'^3 * Y' P_a^2; principal Q'/Y'^3 * Y'^2 P_a^(1/6):
  checks the dyadic bookkeeping sum_{q_k<<P_a/q_C} q_{Ck} ~ P_a^2/q_C and sum_C0 C0^(-1/6) << 1 symbolically.
Run: uv run --no-project --with sympy python3 lemma133_exponents.py
"""
import random, sympy as sp
O, H, S0, B0, za, v, lb, e, A0, N0, l, d = sp.symbols('O H S0 B0 z_a v l_b e A0 N0 l d', real=True)
lp = l - d; Mp = 1 - l - 2*d; Ns = 1 + lp
Td = 2*H + 2*A0 + 2*za - Ns - N0 - 3*B0
DH = Mp - O - H
lhs = O/2 + DH + S0 + B0 + Td/4
rhs = (2*Mp + 2*za - 1 - lp + 2*DH)/4 + (2*A0 - N0 + B0 + 4*S0)/4
ident = sp.simplify(lhs - rhs) == 0
print("identity (l.8355-8359):", ident)
fails = 0 if ident else 1
pos = lambda x: max(x, 0.0)
def E(Ov, Hv, S0v, B0v, zav, vv, lbv, ev, A0v, N0v, lv, dv):
    lpv = lv - dv; Td_ = 2*Hv + 2*A0v + 2*zav - (1 + lpv) - N0v - 3*B0v
    u = min(vv, zav, (vv + zav)/3)
    return Ov/2 + Hv - S0v - B0v + zav - u - lbv - 2*ev/3 - pos(Td_ - vv - 3*lbv - ev)/2
random.seed(7); worst1 = -1e9; worst2 = -1e9; n = 0
for _ in range(400000):
    lv = random.uniform(0, 0.2); dv = random.uniform(0, lv); Mpv = 1 - lv - 2*dv
    Ov = random.uniform(0, Mpv); Hv = random.uniform(0, Mpv - Ov)
    A0v = random.uniform(0, Ov/2); N0v = random.uniform(0, A0v)
    S0v = random.uniform(0, 1); B0v = random.uniform(0, 1)
    zav = random.uniform(0, lv - dv); vv = random.uniform(0, 2); lbv = random.uniform(0, 1); ev = random.uniform(0, 1)
    if random.random() < 0.3: S0v = B0v = 0.0
    if random.random() < 0.3: ev = 0.0
    E0 = E(Ov, Hv, S0v, B0v, zav, vv, lbv, 0.0, A0v, N0v, lv, dv)
    Ee = E(Ov, Hv, S0v, B0v, zav, vv, lbv, ev, A0v, N0v, lv, dv)
    worst1 = max(worst1, Ee - (E0 - ev/6))
    worst2 = max(worst2, E0 - (Mpv + pos(dv + 5*lv - 1)/4))
    n += 1
ok1 = worst1 <= 1e-12; ok2 = worst2 <= 1e-12
print("E(e) <= E(0) - e/6 on %d samples: max excess %.3e %s" % (n, worst1, "ok" if ok1 else "FAIL"))
print("E(0) <= M' + (d+5l-1)_+/4 on %d samples: max excess %.3e %s" % (n, worst2, "ok" if ok2 else "FAIL"))
fails += (not ok1) + (not ok2)
# (13.7): dyadic bookkeeping
Pa, qC, Yp, Qp, C0 = sp.symbols('P_a q_C Y Q C_0', positive=True)
nonpr = sp.simplify((Qp/Yp**3)*Yp*Pa**2 - (Qp/Yp)*(Pa**2/Yp))
princ = sp.simplify((Qp/Yp**3)*Yp**2*Pa**sp.Rational(1, 6) - (Qp/Yp)*Pa**sp.Rational(1, 6))
zero = sp.simplify((Qp/Yp**3)*Qp**0*Yp*Yp - Qp/Yp)
geo = sp.summation(sp.Rational(2)**(-sp.Symbol('k', integer=True)/6), (sp.Symbol('k', integer=True), 0, sp.oo))
print("(13.7) nonprincipal %s, principal %s, zero-frequency %s, sum_dyadic C0^(-1/6) = %s" % (nonpr == 0, princ == 0, zero == 0, sp.nsimplify(geo)))
fails += (nonpr != 0) + (princ != 0) + (zero != 0)
print("RESULT:", "PASS" if fails == 0 else "FAIL (%d)" % fails)
