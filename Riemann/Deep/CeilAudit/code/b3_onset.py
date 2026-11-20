# ceilaudit B3: onset theta_c of E1 below 1/2 if the moment holds for nu < (3 theta - 1)/4 (exact profile, no Taylor step).
# Frozen SL/EXP-010 detector (as nu12.py, read-only copy of its construction), R scanned; h(theta;kappa) > 0 iff kappa > k*(theta).
import json, sys
from flint import fmpq, fmpq_poly, arb, ctx
from mpmath import mp, mpf, cot, sqrt, findroot
ctx.prec = 4000; mp.dps = 30
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
pre = {}
def kappa(nu, Rs):
    best = None
    for R in Rs:
        if R not in pre:
            t = 2*R; a1, b1 = ie(Q2, t); W = R*Qp + dQ; a2, b2 = ie(W*W, t); pre[R] = (a1, b1, a2, b2, t)
        a1, b1, a2, b2, t = pre[R]
        c = arb((CP/nu)*a1 + nu*BP*a2)*arb(t).exp() + arb(fmpq(1,2) + (CP/nu)*b1 + nu*BP*b2)
        kap = 1 - c.log()/arb(R)
        if best is None or kap > best[1]: best = (R, kap)
    return best
cth = lambda t: 2 - t/2 - cot(t/sqrt(2))/sqrt(2)
kstar = lambda t: 1 - 2/(2 - cth(t))
def theta_h(k):   # least theta with kappa > k*(theta): solve k*(theta) = k
    return findroot(lambda t: kstar(t) - k, mpf('0.495'))
Rs = [fmpq(780 + 5*i, 100) for i in range(0, 25)]   # 7.80 .. 9.00
_c = kappa(fmpq(3, 25), Rs); print("check: nu=3/25 -> R", _c[0], "kappa", _c[1].str(15, radius=False), "(E1: R=177/20, kappa=0.0844779079697)", flush=True)
rows = []
for num in range(1180, 1251, 5):
    nu = fmpq(num, 10000); R, k = kappa(nu, Rs)
    kf = mpf(k.str(25, radius=False)); th_m = (4*mpf(num)/10000 + 1)/3; th_d = theta_h(kf)
    rows.append((num, R, kf, th_m, th_d, max(th_m, th_d)))
    print(f"nu={num}/10000 R={R} kappa={float(kf):.10f} theta_moment>(4nu+1)/3={float(th_m):.6f} theta_detector={float(th_d):.6f} -> onset {float(max(th_m, th_d)):.6f}", flush=True)
b = min(rows, key=lambda x: x[5])
print(f"BEST onset on the nu grid: theta_c = {float(b[5]):.6f} at nu = {b[0]}/10000 (R = {b[1]}); h(1/2; kappa) there would be positive for every theta > theta_c")
