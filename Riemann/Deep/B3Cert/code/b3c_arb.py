# b3cert step (5): Arb certificates. kappa(121/1000) at R = 44/5 with the frozen SL/EXP-010 detector (construction copied from
# Riemann/Deep/SantiAudit/nu12.py, read-only), Wang's c(theta), k*(theta) = 1 - 2/(2 - c(theta)), h(theta;k) = (3 + k - sqrt((1-k)(9 - k - 8c)))/4.
import json, sys
from flint import fmpq, fmpq_poly, arb, ctx
ctx.prec = 4000
D = json.load(open(sys.argv[1]))
ps = [p for p in D['parameter_sets'] if p['nu'] == '999/10000'][0]
q = lambda s: fmpq(*map(int, s.split('/'))) if '/' in s else fmpq(int(s))
xs = [q(x) for x in ps['x_j']]; r = q(ps['r_P'])
u = fmpq_poly([1, -2]); T = [fmpq_poly([1]), u]
for n in range(1, 2*len(xs)+1): T.append(2*u*T[n] - T[n-1])
Qp = fmpq_poly([1, -1])
for j, x in enumerate(xs, start=1): Qp += x*(T[2*j+1] - u)
S = fmpq_poly([0])
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
def kappa(nu, R):
    t = 2*R; a1, b1 = ie(Q2, t); W = R*Qp + dQ; a2, b2 = ie(W*W, t)
    c = arb((CP/nu)*a1 + nu*BP*a2)*arb(t).exp() + arb(fmpq(1,2) + (CP/nu)*b1 + nu*BP*b2)
    return 1 - c.log()/arb(R), c
# admissibility of the detector (exact rationals): P(0)=0, P(1)=1, Q(0)=1, Q(y)+Q(1-y) = 1 identically (Q = integral form as SL: check via Qp)
print("P(0) =", P(fmpq(0)), " P(1) =", P(fmpq(1)), " Qp(0) =", Qp(fmpq(0)), " Qp(y)+Qp(1-y)-1 == 0:", (Qp + Qp(fmpq_poly([1, -1])) - 1) == 0, " deg Qp =", Qp.degree(), flush=True)
k3, c3 = kappa(fmpq(3, 25), fmpq(177, 20)); print("check E1: nu=3/25 R=177/20 kappa =", k3.str(30), flush=True)
nu = fmpq(121, 1000); R = fmpq(44, 5); k, c = kappa(nu, R)
print("nu = 121/1000, R = 44/5: c =", c.str(30), " kappa =", k.str(30), flush=True)
s2 = arb(2).sqrt()
cth = lambda th: 2 - th/2 - (th/s2).cot()/s2
kst = lambda th: 1 - 2/(2 - cth(th))
h = lambda th, kk: (3 + kk - ((1 - kk)*(9 - kk - 8*cth(th))).sqrt())/4
th0 = arb(fmpq(371, 750))
print("moment edge: (3*371/750 - 1)/4 == 121/1000:", (3*fmpq(371, 750) - 1)/4 == nu, flush=True)
for name, th in [("371/750", th0), ("4947/10000", arb(fmpq(4947, 10000))), ("1/2", arb(fmpq(1, 2)))]:
    print(f"theta = {name}: c = {cth(th).str(25)}, k* = {kst(th).str(25)}, kappa - k* = {(k - kst(th)).str(20)} (>0: {(k - kst(th)) > 0}), "
          f"h = {h(th, k).str(25)} (>0: {h(th, k) > 0})", flush=True)
# detector threshold: bisection for k*(theta) = kappa (k* strictly decreasing since c' = cot^2(theta/sqrt2)/2 > 0, Wang l.129)
lo, hi = fmpq(49, 100), fmpq(371, 750)
assert kst(arb(lo)) > k and kst(arb(hi)) < k
for _ in range(60):
    mid = (lo + hi)/2
    if kst(arb(mid)) > k: lo = mid
    else: hi = mid
print("theta_det in [", arb(lo).str(15), ",", arb(hi).str(15), "]  < 371/750 = 0.494666...:", hi < fmpq(371, 750), flush=True)
