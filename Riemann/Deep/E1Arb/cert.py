# e1arb: Arb certificate of h(theta; kappa) for one detector in E1's class. Exact rationals (flint fmpq, fmpq_poly) for P, Q, nu, R, theta,
# C_P, B_P and the polynomial parts of int_0^1 e^{2Rv}F(v)dv = e^{2R}G1 + G0 (repeated integration by parts, finite, exact); arb balls only
# for e^{2R}, log, sqrt, cot. Every admissibility constraint of E1 (REPORT.md section 1) is checked exactly.
# usage: cert.py par.json [prec]   par: {"theta","nu","R","r","N","x":[...]} rationals as "a/b" strings
import json, sys
from flint import fmpq, fmpq_poly, arb, ctx
def q(s):
    s = str(s); a, b = (s.split('/') + ['1'])[:2]; return fmpq(int(a), int(b))
def buildQ(xs):
    u = fmpq_poly([1, -2]); T = [fmpq_poly([1]), u]
    for k in range(1, 2*len(xs) + 1): T.append(2*u*T[k] - T[k-1])
    Q = fmpq_poly([1, -1])
    for j, x in enumerate(xs, start=1): Q += x*(T[2*j+1] - u)
    return Q
def buildP(r, N):
    S = fmpq_poly([0]); f = fmpq(1)
    for k in range(1, 2*N):
        f = f/k
        if k % 2: S += fmpq_poly([0]*k + [f])
    return S(fmpq_poly([0, r]))*(1/S(r))
def checks(th, nu, R, P, Q, uniform=False):  # uniform: theta = 1/2, nu < 1/8 = (3/2-1)/4 is admissible at every theta in (1/2,1), h(1/2;kappa) <= h(theta;kappa)
    return {'0<nu': nu > 0, 'nu<(3theta-1)/4': nu < (3*th - 1)/4, 'nu<1/2': nu < fmpq(1, 2), '1/2<theta<1': fmpq(1, 2) < th < 1 or (uniform and th == fmpq(1, 2)), 'R>0': R > 0,
            'Q(0)=1': Q(fmpq(0)) == 1, 'Q(1)=0': Q(fmpq(1)) == 0, 'Q(y)+Q(1-y)=1': Q + Q(fmpq_poly([1, -1])) == fmpq_poly([1]),
            'P(0)=0': P(fmpq(0)) == 0, 'P(1)=1': P(fmpq(1)) == 1}
I01 = lambda F: F.integral()(fmpq(1)) - F.integral()(fmpq(0))
def ie(F, t):  # int_0^1 e^{tv}F(v)dv = e^t*G1 + G0, exact rationals
    G1 = fmpq(0); G0 = fmpq(0); Fj = F; j = 0
    while Fj.degree() >= 0:
        s = fmpq((-1)**j)/t**(j+1); G1 += s*Fj(fmpq(1)); G0 -= s*Fj(fmpq(0)); Fj = Fj.derivative(); j += 1
    return G1, G0
def kappa(nu, R, P, Q):
    dP = P.derivative(); CP, BP = I01(dP*dP), I01(P*P); t = 2*R
    a1, b1 = ie(Q*Q, t); W = R*Q + Q.derivative(); a2, b2 = ie(W*W, t); et = arb(t).exp()
    c = arb((CP/nu)*a1 + nu*BP*a2)*et + arb(fmpq(1, 2) + (CP/nu)*b1 + nu*BP*b2)
    return c, 1 - c.log()/arb(R), {'CP': CP, 'BP': BP, 'W2': arb(a1)*et + arb(b1), 'W2p': arb(a2)*et + arb(b2)}
def cwang(th):
    s2 = arb(2).sqrt(); x = arb(th)/s2; return 2 - arb(th)/2 - (x.cos()/x.sin())/s2
def hball(th, k):
    ct = cwang(th); rad = (1 - k)*(9 - k - 8*ct); assert rad > 0 and k < 1 and ct < 1
    return (3 + k - rad.sqrt())/4
def floorcert(b, d=10):  # largest multiple of 10^-d certified below the ball b
    v = fmpq(int((b.mid()*10**d).floor().unique_fmpz()) - 1, 10**d)
    while b > arb(v + fmpq(1, 10**d)): v += fmpq(1, 10**d)
    assert b > arb(v); return v
def load(fn):
    D = json.load(open(fn)); th, nu, R, r = q(D['theta']), q(D['nu']), q(D['R']), q(D['r'])
    return th, nu, R, buildP(r, int(D['N'])), buildQ([q(x) for x in D['x']]), D
if __name__ == "__main__":
    ctx.prec = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
    th, nu, R, P, Q, D = load(sys.argv[1]); chk = checks(th, nu, R, P, Q, D.get('uniform', False)); assert all(chk.values()), chk
    c, k, X = kappa(nu, R, P, Q); h = hball(th, k); ct = cwang(th); lb3 = (ct + 2*k)/3; hb, kb = floorcert(h), floorcert(k)
    print("checks:", ' '.join(kk for kk, vv in chk.items() if vv), f"| deg Q {Q.degree()}, deg P {P.degree()}, prec {ctx.prec}")
    print(f"theta {th} nu {nu} R {R} r {D['r']} N {D['N']} m {len(D['x'])}")
    print("C_P =", arb(X['CP']).str(20), " B_P =", arb(X['BP']).str(20), " W2 =", X['W2'].str(20), " W2' =", X['W2p'].str(20))
    print("c(P,Q,R,nu) =", c.str(30)); print("kappa =", k.str(30)); print("c(theta) Wang =", ct.str(25))
    print("h(theta;kappa) =", h.str(30)); print("(c+2kappa)/3 =", lb3.str(25))
    print(f"CERTIFIED: kappa > {kb} = {float(kb):.10f};  h(theta;kappa) > {hb} = {float(hb):.10f}")
    json.dump({'theta': str(th), 'nu': str(nu), 'kappa_ball': k.str(40), 'h_ball': h.str(40), 'kappa_gt': str(kb), 'h_gt': str(hb),
               'c_theta': ct.str(30), 'lb3_ball': lb3.str(30), 'prec': ctx.prec}, open(sys.argv[1].replace('.json', '.cert.json'), 'w'), indent=1)
