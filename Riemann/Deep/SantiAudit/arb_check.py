# Recompute c(P,Q,R,nu), kappa, h(theta;kappa) for the paper's frozen detector, exact rationals + Arb.
import json, sys
from flint import fmpq, fmpq_poly, arb, ctx
ctx.prec = 4000
D = json.load(open(sys.argv[1]))
pick = sys.argv[2] if len(sys.argv) > 2 else "499/10000"
theta = fmpq(527, 1000) if len(sys.argv) <= 3 else fmpq(*map(int, sys.argv[3].split('/')))

def q(s):
    s = str(s)
    if '/' in s:
        a, b = s.split('/'); return fmpq(int(a), int(b))
    return fmpq(int(s))

def integ_exp(F, t):
    # exact: int_0^1 e^{t v} F(v) dv = e^t * G(1) - G(0), G = sum_j (-1)^j F^(j) / t^(j+1)
    G1 = fmpq(0); G0 = fmpq(0); Fj = F; j = 0
    while Fj.degree() >= 0:
        s = fmpq((-1) ** j) / t ** (j + 1)
        G1 += s * Fj(fmpq(1)); G0 += s * Fj(fmpq(0))
        Fj = Fj.derivative(); j += 1
    return G1, -G0  # coefficient of e^t, constant

def integ01(F):
    A = F.integral()
    return A(fmpq(1)) - A(fmpq(0))

for ps in D['parameter_sets']:
    if ps['nu'] != pick: continue
    nu, R, r = q(ps['nu']), q(ps['R']), q(ps['r_P'])
    xs = [q(x) for x in ps['x_j']]
    u = fmpq_poly([1, -2])
    T = [fmpq_poly([1]), u]
    for n in range(1, 2 * len(xs) + 1):
        T.append(2 * u * T[n] - T[n - 1])
    Qp = fmpq_poly([1, -1])
    for j, x in enumerate(xs, start=1):
        Qp += x * (T[2 * j + 1] - u)
    one_minus_y = fmpq_poly([1, -1])
    sym = Qp(one_minus_y) + Qp
    print("deg Q", Qp.degree(), "Q(0)", Qp(fmpq(0)), "Q(1)", Qp(fmpq(1)), "Q(y)+Q(1-y) == 1:", sym == fmpq_poly([1]))
    fact = 1; S = fmpq_poly([0])
    for i in range(7):
        k = 2 * i + 1
        fact = 1
        for m in range(2, k + 1): fact *= m
        S += fmpq_poly([0] * k + [fmpq(1, fact)])
    Sr = S(fmpq_poly([0, r]))
    P = Sr / S(r) if False else Sr * (fmpq(1) / S(r))
    print("deg P", P.degree(), "P(0)", P(fmpq(0)), "P(1)", P(fmpq(1)))
    dP = P.derivative()
    CP = integ01(dP * dP); BP = integ01(P * P); PP = integ01(P * dP)
    print("int P P' =", PP)
    t = 2 * R
    a1, b1 = integ_exp(Qp * Qp, t)
    W = R * Qp + Qp.derivative()
    a2, b2 = integ_exp(W * W, t)
    A = (CP / nu) * a1 + nu * BP * a2
    B = fmpq(1, 2) + (CP / nu) * b1 + nu * BP * b2
    c = arb(A) * arb(t).exp() + arb(B)
    kappa = 1 - c.log() / arb(R)
    # direct (unseparated) check of int w w' = -1/2 : e^{2R}Q(1)^2/2 - Q(0)^2/2
    print("nu", nu, "R", R, "r", r)
    print("c =", c.str(30))
    print("kappa =", kappa.str(40))
    kp = arb(fmpq(71589990975836346694543352023356087639, 2 * 10 ** 39))
    print("paper kappa lower 7158999.../(2e39) =", kp.str(40), " kappa - paper_lb > 0:", (kappa - kp) > 0)
    th = arb(theta); s2 = arb(2).sqrt()
    x = th / s2
    cth = 2 - th / 2 - (x.cos() / x.sin()) / s2
    h = (3 + kappa - ((1 - kappa) * (9 - kappa - 8 * cth)).sqrt()) / 4
    print("theta", theta, "c(theta) =", cth.str(30))
    print("onset k* = 1-2/(2-c) =", (1 - 2 / (2 - cth)).str(30))
    print("h(theta;kappa) =", h.str(30), " h > 0.0005947542001:", h > arb(fmpq(5947542001, 10 ** 13)))
    print("(c+2k)/3 =", ((cth + 2 * kappa) / 3).str(20))
    print("nuA(theta) = min(1/2, 2th-1, (5th-2)/6) =", min(fmpq(1, 2), 2 * theta - 1, (5 * theta - 2) / 6), " nu < nuA:", nu < min(fmpq(1, 2), 2 * theta - 1, (5 * theta - 2) / 6))
