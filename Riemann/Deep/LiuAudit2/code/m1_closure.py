"""liuaudit2 T11: Lemma 10.6 / §10.4 exact identities (number_fields_v4.tex l.4903-4910) from the child centers (l.4612-4618),
F = N + V (l.3930), r = N - 3l (l.4026); then the closure inequalities (l.4918-4924) from the stated allowances:
M' <= Mhat + 14 eta (row allowance), Q' <= Q + lhat + that + ghat - thetahat + eta, F' >= Fhat, z' <= z0,
ghat - thetahat >= -eta (l.4899), Ai, t, j, V >= 0, |lhat - l| <= eta, Mhat <= M - 2(V + lhat) + 3 eta (l.4927).
Run: uv run --no-project --with sympy python3 m1_closure.py
"""
import sympy as sp, random
M, Q, V, l, Ai, B, t, g, v, th, j, N = sp.symbols('M Q V l A_i B t g v theta j N', real=True)
r = N - 3*l; F = N + V
Nc = r - Ai - B - t - g - v; Vc = B + th + v + j; Fc = Nc + Vc
Mc = M - 4*l - 2*Ai - 2*t - 2*g + 2*th - 2*V + j
id1 = sp.simplify(Fc - Mc - (Q + l + t + g - th) - (F - M - Q + Ai + V)) == 0
id2 = sp.simplify(4*Fc - 3*Mc - (4*F - 3*M + 2*Ai + 2*t + 2*(g - th) + 2*V + j)) == 0
print("identity 1 (l.4904-4906):", id1, "  identity 2 (l.4907-4908):", id2)
# inequalities: worst case over the allowances (linear; random extreme-biased sampling, plain floats)
fFc = sp.lambdify((M, Q, V, l, Ai, B, t, g, v, th, j, N), Fc)
fMc = sp.lambdify((M, Q, V, l, Ai, B, t, g, v, th, j, N), Mc)
random.seed(3); worst = [-1e9, -1e9]
for _ in range(200000):
    eta = random.uniform(0, 0.01)
    Mv, Qv, Vv, lv = random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 0.2)
    Aiv, Bv, tv, vv, jv, Nv = (random.uniform(0, 1) for _ in range(6))
    thv = random.uniform(0, 1); gv = thv - eta*random.random() + random.uniform(0, 1)*random.choice([0, 1])
    Fh = fFc(Mv, Qv, Vv, lv, Aiv, Bv, tv, gv, vv, thv, jv, Nv); Mh = fMc(Mv, Qv, Vv, lv, Aiv, Bv, tv, gv, vv, thv, jv, Nv)
    Fp = Fh + random.uniform(0, 1)*random.choice([0, 1]); Mp = Mh + 14*eta*random.random()
    Qp = Qv + lv + tv + gv - thv + eta*random.random()
    z0 = random.uniform(0, 1); zp = z0*random.random()
    Fv = Nv + Vv
    worst[0] = max(worst[0], (Fv - Mv - Qv - z0 - 15*eta) - (Fp - Mp - Qp - zp))
    worst[1] = max(worst[1], (4*Fv - 3*Mv - 6*z0 - 44*eta) - (4*Fp - 3*Mp - 6*zp))
ok = worst[0] <= 1e-12 and worst[1] <= 1e-12
print("closure inequalities (l.4918-4919) under the allowances: max violation %.2e, %.2e -> %s" % (worst[0], worst[1], "ok" if ok else "FAIL"))
print("RESULT:", "PASS" if (id1 and id2 and ok) else "FAIL")
