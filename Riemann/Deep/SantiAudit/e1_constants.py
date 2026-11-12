# E1 constants in Arb: detector (P,Q) of the frozen nu=999/10000 EXP-010 set at nu=3/25, R=177/20 (and nu=1249/10000, R=173/20).
import json
from flint import fmpq, fmpq_poly, arb, ctx
ctx.prec = 4000
D = json.load(open('params_exp010.json'))
ps = [p for p in D['parameter_sets'] if p['nu'] == '999/10000'][0]
q = lambda s: fmpq(*map(int, s.split('/'))) if '/' in s else fmpq(int(s))
xs = [q(x) for x in ps['x_j']]; r = q(ps['r_P'])
u = fmpq_poly([1, -2]); T = [fmpq_poly([1]), u]
for n in range(1, 2*len(xs) + 1): T.append(2*u*T[n] - T[n-1])
Qp = fmpq_poly([1, -1])
for j, x in enumerate(xs, start=1): Qp += x*(T[2*j+1] - u)
S = fmpq_poly([0])
for i in range(7):
    kk = 2*i + 1; f = 1
    for m in range(2, kk + 1): f *= m
    S += fmpq_poly([0]*kk + [fmpq(1, f)])
P = S(fmpq_poly([0, r])) * (fmpq(1)/S(r)); dP = P.derivative()
print("Q: deg", Qp.degree(), "Q(0)=", Qp(fmpq(0)), "Q(1)=", Qp(fmpq(1)), "Q(y)+Q(1-y)==1:", Qp(fmpq_poly([1, -1])) + Qp == fmpq_poly([1]),
      "| P: deg", P.degree(), "P(0)=", P(fmpq(0)), "P(1)=", P(fmpq(1)))
I01 = lambda F: F.integral()(fmpq(1)) - F.integral()(fmpq(0))
CP, BP = I01(dP*dP), I01(P*P)
def ie(F, t):
    G1 = fmpq(0); G0 = fmpq(0); Fj = F; j = 0
    while Fj.degree() >= 0:
        s = fmpq((-1)**j)/t**(j+1); G1 += s*Fj(fmpq(1)); G0 += s*Fj(fmpq(0)); Fj = Fj.derivative(); j += 1
    return G1, -G0
s2 = arb(2).sqrt()
def cth(th):
    x = th/s2; return 2 - th/2 - (x.cos()/x.sin())/s2
def h(th, k): return (3 + k - ((1 - k)*(9 - k - 8*cth(th))).sqrt())/4
for nu, R in [(fmpq(3, 25), fmpq(177, 20)), (fmpq(1249, 10000), fmpq(173, 20))]:
    t = 2*R; a1, b1 = ie(Qp*Qp, t); W = R*Qp + Qp.derivative(); a2, b2 = ie(W*W, t)
    c = arb((CP/nu)*a1 + nu*BP*a2)*arb(t).exp() + arb(fmpq(1, 2) + (CP/nu)*b1 + nu*BP*b2)
    k = 1 - c.log()/arb(R)
    half = arb(fmpq(1, 2)); c5 = cth(half); kstar = 1 - 2/(2 - c5)
    h5 = h(half, k)
    x = half/s2; cp = (x.cos()/x.sin())**2/2
    dh = ((1 - k)/(9 - k - 8*c5)).sqrt()*cp
    print("nu", nu, "R", R, "theta-threshold (4nu+1)/3 =", (4*nu + 1)/3)
    print("  c =", c.str(25), " kappa =", k.str(30))
    print("  c(1/2) =", c5.str(25), " k*(1/2) =", kstar.str(25), " kappa - k*(1/2) =", (k - kstar).str(15))
    print("  h(1/2;kappa) =", h5.str(25), " (lower bound for every theta in (1/2,1): h increasing in theta)")
    print("  dh/dtheta at 1/2 =", dh.str(15), "; c'(1/2) =", cp.str(15))
    for th in ['1/2', '501/1000', '51/100', '52/100', '55/100', '6/10', '3/4']:
        print("   theta", th, " h =", h(arb(q(th)), k).str(12))
