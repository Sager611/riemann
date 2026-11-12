# Re-use the frozen nu=999/10000 detector polynomials (P,Q) at a longer mollifier nu, scanning rational R.
import json, sys
from flint import fmpq, fmpq_poly, arb, ctx
ctx.prec = 4000
D = json.load(open(sys.argv[1])); nu = fmpq(*map(int, sys.argv[2].split('/')))
ps = [p for p in D['parameter_sets'] if p['nu'] == '999/10000'][0]
q = lambda s: fmpq(*map(int, s.split('/'))) if '/' in s else fmpq(int(s))
xs = [q(x) for x in ps['x_j']]; r = q(ps['r_P'])
u = fmpq_poly([1, -2]); T = [fmpq_poly([1]), u]
for n in range(1, 2*len(xs)+1): T.append(2*u*T[n] - T[n-1])
Qp = fmpq_poly([1, -1])
for j, x in enumerate(xs, start=1): Qp += x*(T[2*j+1] - u)
S = fmpq_poly([0]); f = 1
for i in range(7):
    kk = 2*i+1; f = 1
    for m in range(2, kk+1): f *= m
    S += fmpq_poly([0]*kk + [fmpq(1, f)])
P = S(fmpq_poly([0, r])) * (fmpq(1)/S(r)); dP = P.derivative()
I01 = lambda F: F.integral()(fmpq(1)) - F.integral()(fmpq(0))
CP, BP = I01(dP*dP), I01(P*P)
def ie(F, t):
    G1 = fmpq(0); G0 = fmpq(0); Fj = F; j = 0
    while Fj.degree() >= 0:
        s = fmpq((-1)**j)/t**(j+1); G1 += s*Fj(fmpq(1)); G0 += s*Fj(fmpq(0)); Fj = Fj.derivative(); j += 1
    return G1, -G0
Q2 = Qp*Qp; dQ = Qp.derivative()
best = None
for i in range(0, 81):
    R = fmpq(700 + 5*i, 100)   # 7.00 .. 11.00
    t = 2*R; a1, b1 = ie(Q2, t); W = R*Qp + dQ; a2, b2 = ie(W*W, t)
    c = arb((CP/nu)*a1 + nu*BP*a2)*arb(t).exp() + arb(fmpq(1,2) + (CP/nu)*b1 + nu*BP*b2)
    kap = 1 - c.log()/arb(R)
    if best is None or kap > best[1]: best = (R, kap, c)
print("nu", nu, "best R", best[0], "kappa", best[1].str(25), "c", best[2].str(15))
k5 = arb(fmpq(7662355506857570169601771519720957, 10**35))   # k*(1/2) = 1 - 2/(2 - c(1/2)) (rounded up below)
print("k*(1/2) ~ 0.0766235550685757017; kappa > k*(1/2):", best[1] > arb(fmpq(766235550685758, 10**16)))
print("nu < (3*theta-1)/4 for theta>1/2 needs nu <= 1/8:", nu <= fmpq(1, 8))
