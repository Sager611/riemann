# e1rate: kappa at the optimized (R, Q) recomputed in mpmath (P = sinh(l x)/sinh l and its degree-15 odd Taylor polynomial), then the h(theta) table.
import json, mpmath as mp
mp.mp.dps = 30
D = {float(k): v for k, v in json.load(open('kappa_opt.json')).items()}
def Q_and_dQ(xs, y):
    u = 1-2*y; q = (1+u)/2; dq = mp.mpf(-1)
    for j, x in enumerate(xs, start=1):
        n = 2*j+1; T = mp.chebyt(n, u); dT = n*mp.chebyu(n-1, u)
        q += x*(T-u); dq += x*(dT-1)*(-2)
    return q, dq
def kap(nu, R, xs, poly=False):
    nu = mp.mpf(nu); R = mp.mpf(R); xs = [mp.mpf(x) for x in xs]
    W2 = mp.quad(lambda v: mp.e**(2*R*v)*Q_and_dQ(xs, v)[0]**2, [0, 1])
    W2p = mp.quad(lambda v: mp.e**(2*R*v)*(R*Q_and_dQ(xs, v)[0]+Q_and_dQ(xs, v)[1])**2, [0, 1])
    a = W2/nu; b = nu*W2p; l = mp.sqrt(b/a)
    if not poly:
        c = mp.mpf(1)/2 + mp.sqrt(a*b)*mp.coth(l)
    else:  # P = S(l x)/S(l), S = odd Taylor polynomial of sinh to degree 15 (an admissible polynomial: P(0)=0, P(1)=1)
        S = lambda t: sum(t**(2*i+1)/mp.factorial(2*i+1) for i in range(8))
        dS = lambda t: sum(t**(2*i)/mp.factorial(2*i) for i in range(8))
        CP = mp.quad(lambda x: (l*dS(l*x)/S(l))**2, [0, 1]); BP = mp.quad(lambda x: (S(l*x)/S(l))**2, [0, 1])
        c = mp.mpf(1)/2 + a*CP + b*BP
    return 1-mp.log(c)/R
s2 = mp.sqrt(2)
def cth(th): x = th/s2; return 2-th/2-mp.cot(x)/s2
def h(th, k): return (3+k-mp.sqrt((1-k)*(9-k-8*cth(th))))/4
def LB(th, k): c = cth(th); return max(mp.mpf(0), c, (c+2*k)/3, h(th, k))
K = {}
for nu, (k, p) in sorted(D.items()):
    ke = kap(nu, p[0], p[1:]); kp = kap(nu, p[0], p[1:], poly=True); K[round(nu, 6)] = ke
    print(f"nu={nu:.6f} kappa_float={k:.10f} kappa_mp(sinh P)={mp.nstr(ke, 12)} kappa_mp(deg-15 poly P)={mp.nstr(kp, 12)} R={p[0]:.5f}")
def kof(nu): return K[round(nu, 6)]
print("theta | nu_E1 | nu_B(SL) | old nu = max | LB_old | nu_K (new) | LB_new | gain")
rows = [(0.51, None), (0.6, None), (0.75, None), (0.82, None), (0.85, None), (0.9, None), (0.95, None), (0.99, None)]
for th, _ in rows:
    nE = (3*th-1)/4; nB = min(0.5, 17/33*(2*th-1)); old = max(nE, nB)
    nK = min(0.5, max(nE, 4/7*(2*th-1)))
    ko = kof(old) if round(old, 6) in K else None; kn = kof(nK) if round(nK, 6) in K else None
    lo = LB(mp.mpf(th), ko) if ko is not None else None; ln = LB(mp.mpf(th), kn) if kn is not None else None
    print(f"{th} | {nE:.6f} | {nB:.6f} | {old:.6f} | {mp.nstr(lo, 9) if lo is not None else '-'} | {nK:.6f} | {mp.nstr(ln, 9) if ln is not None else '-'} | {mp.nstr(ln-lo, 6) if (lo is not None and ln is not None) else '-'}")
th = mp.mpf(0.99); print("uncapped (conditional) theta=0.99, nu=0.56: kappa", mp.nstr(kof(0.56), 10), "LB", mp.nstr(LB(th, kof(0.56)), 9))
for th in [0.51, 0.6, 0.75, 0.82, 0.85, 0.9, 0.95, 0.99]: print("c(theta)", th, mp.nstr(cth(mp.mpf(th)), 10))
