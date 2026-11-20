# b3audit: independent recomputation of kappa(nu, R) = 1 - log(c)/R, c = Q(0)^2 + nu^{-1} int_0^1 int_0^1 (w_R P' + nu w_R' P)^2 du dv,
# w_R(v) = e^{Rv} Q(v) (SL l. 174), for the frozen EXP-010 detector (nu = 999/10000 set of params_exp010.json), and of Wang's c(theta),
# k*(theta) = 1 - 2/(2 - c(theta)), h(theta; k). Method differs from b3c_arb.py (exact rational antiderivatives in Arb): floating mpmath at
# 50 digits, Q and Q' from the trigonometric form T_n(cos phi) = cos(n phi), T_n' = n sin(n phi)/sin(phi) (no polynomial coefficients),
# Gauss-Legendre quadrature on 24 subintervals; the cross term int w_R w_R' and int P P' are integrated numerically (Q(1) = 0 not used).
import json, sys, mpmath as mp
mp.mp.dps = 50
D = json.load(open(sys.argv[1]))
ps = [p for p in D['parameter_sets'] if p['nu'] == '999/10000'][0]
def q(s):
    s = str(s)
    if '/' in s: a, b = s.split('/'); return mp.mpf(int(a))/int(b)
    return mp.mpf(s)
xs = [q(x) for x in ps['x_j']]; rP = q(ps['r_P'])
S = lambda t: mp.fsum(t**(2*i+1)/mp.factorial(2*i+1) for i in range(7))
dS = lambda t: mp.fsum(t**(2*i)/mp.factorial(2*i) for i in range(7))
P = lambda x: S(rP*x)/S(rP); dP = lambda x: rP*dS(rP*x)/S(rP)
def QdQ(x):   # e^{i n phi} by the recurrence z^{n+2} = z^n z^2 (stable, |z| = 1)
    u = 1 - 2*x; ph = mp.acos(u); s = mp.sin(ph); Qv = 1 - x; dQv = mp.mpf(-1)
    z = mp.expj(ph); z2 = z*z; zn = z*z2
    for j, c in enumerate(xs, start=1):
        n = 2*j + 1; Qv += c*(zn.real - u); dQv += -2*c*(n*zn.imag/s - 1); zn *= z2
    return Qv, dQv
from mpmath.calculus.quadrature import GaussLegendre
GL = GaussLegendre(mp.mp); DEG = int(sys.argv[2]) if len(sys.argv) > 2 else 7
NODES = []
for i in range(24):
    NODES += GL.transform_nodes(GL.calc_nodes(DEG, mp.mp.prec), mp.mpf(i)/24, mp.mpf(i+1)/24)
print("Gauss-Legendre nodes:", len(NODES), "degree", DEG, flush=True)
gl = lambda f: mp.fsum(w*f(x) for x, w in NODES)
def glv(f):
    acc = None
    for x, w in NODES:
        v = f(x); acc = [w*t for t in v] if acc is None else [a + w*t for a, t in zip(acc, v)]
    return acc
CP, BP, PdP = gl(lambda x: dP(x)**2), gl(lambda x: P(x)**2), gl(lambda x: P(x)*dP(x))
print("P(0) =", mp.nstr(P(0), 5), " P(1) =", mp.nstr(P(1), 30), " int P P' =", mp.nstr(PdP, 30), flush=True)
Q0, _ = QdQ(mp.mpf('1e-30')); Q1, _ = QdQ(1 - mp.mpf('1e-30')); Qm, _ = QdQ(mp.mpf('0.3')); Qm2, _ = QdQ(mp.mpf('0.7'))
print("Q(0+) =", mp.nstr(Q0, 20), " Q(1-) =", mp.nstr(Q1, 20), " Q(.3)+Q(.7)-1 =", mp.nstr(Qm + Qm2 - 1, 5), flush=True)
def kappa(nu, R):
    def ints(v):
        Qv, dQv = QdQ(v); e = mp.exp(2*R*v); W = R*Qv + dQv
        return Qv*Qv*e, Qv*W*e, W*W*e
    IQQ, IQW, IWW = glv(ints)
    c = 1 + (CP*IQQ + 2*nu*PdP*IQW + nu*nu*BP*IWW)/nu
    return 1 - mp.log(c)/R, c, IQW
k3, c3, _ = kappa(mp.mpf(3)/25, mp.mpf(177)/20); print("check E1: nu=3/25 R=177/20 kappa =", mp.nstr(k3, 30), flush=True)
k, c, IQW = kappa(mp.mpf(121)/1000, mp.mpf(44)/5)
print("nu=121/1000 R=44/5: c =", mp.nstr(c, 30), " kappa =", mp.nstr(k, 30), " int w_R w_R' =", mp.nstr(IQW, 25), flush=True)
cth = lambda th: 2 - th/2 - mp.cot(th/mp.sqrt(2))/mp.sqrt(2)
kst = lambda th: 1 - 2/(2 - cth(th))
h = lambda th, kk: (3 + kk - mp.sqrt((1 - kk)*(9 - kk - 8*cth(th))))/4
th0 = mp.mpf(371)/750
print("theta0=371/750: c =", mp.nstr(cth(th0), 25), " k* =", mp.nstr(kst(th0), 25), " kappa-k* =", mp.nstr(k - kst(th0), 20),
      " h =", mp.nstr(h(th0, k), 25), flush=True)
thd = mp.findroot(lambda t: kst(t) - k, mp.mpf('0.4946'))
print("theta_det (k*(theta) = kappa, findroot) =", mp.nstr(thd, 20), " < 371/750:", thd < th0, flush=True)
