# E1Audit: independent recomputation of c(P,Q,R,nu), kappa1, k*(1/2), h(1/2;kappa1) for the EXP-010 nu=999/10000 detector at nu=3/25, R=177/20.
# Route (differs from E1's exact-rational monomial route): Q evaluated in Chebyshev form T_n(t)=cos(n acos t), T_n'(t)=n sin(n th)/sin th,
# Gauss-Legendre quadrature in mpmath at two node counts and two precisions; downstream c(theta), k*, h in Arb (python-flint) with the
# quadrature discrepancy (times 10) as the radius of c.
import json, sys
from mpmath import mp, mpf, cos, sin, acos, exp, log, sqrt, cot
D = json.load(open('../SantiAudit/params_exp010.json'))
ps = [p for p in D['parameter_sets'] if p['nu'] == '999/10000'][0]
fr = lambda s: mpf(int(s.split('/')[0])) / int(s.split('/')[1]) if '/' in s else mpf(int(s))
def comp(dps, N, nu_s, R_s):
    mp.dps = dps
    xs = [fr(x) for x in ps['x_j']]; rP = fr(ps['r_P']); nu = fr(nu_s); R = fr(R_s)
    S = lambda u: sum(u ** (2 * i + 1) / mp.factorial(2 * i + 1) for i in range(7))
    dS = lambda u: sum(u ** (2 * i) / mp.factorial(2 * i) for i in range(7))
    P = lambda u: S(rP * u) / S(rP); dP = lambda u: rP * dS(rP * u) / S(rP)
    def QdQ(v):
        t = 1 - 2 * v; th = acos(t)
        Q = 1 - v + sum(x * (cos((2 * j + 1) * th) - t) for j, x in enumerate(xs, 1))
        dQdt = sum(x * ((2 * j + 1) * sin((2 * j + 1) * th) / sin(th) - 1) for j, x in enumerate(xs, 1))
        return Q, -1 - 2 * dQdt
    nodes = mp.gauss_legendre(N) if hasattr(mp, 'gauss_legendre') else None
    from mpmath.calculus.quadrature import GaussLegendre
    gl = GaussLegendre(mp); deg = max(1, (N - 1).bit_length())
    pts = gl.calc_nodes(deg, mp.prec)  # nodes on [-1,1] with weights; 3*2^(deg-1) points
    I1 = I2 = CP = BP = mpf(0)
    for (x, wgt) in pts:
        v = (x + 1) / 2; ww = wgt / 2
        Q, dQ = QdQ(v); e2 = exp(2 * R * v)
        I1 += ww * e2 * Q * Q; I2 += ww * e2 * (R * Q + dQ) ** 2
        CP += ww * dP(v) ** 2; BP += ww * P(v) ** 2
    c = mpf(1) / 2 + CP * I1 / nu + nu * BP * I2
    Q0, _ = QdQ(mpf('1e-30')); Q1, _ = QdQ(1 - mpf('1e-30')); Qh, _ = QdQ(mpf(1) / 3); Qh2, _ = QdQ(mpf(2) / 3)
    return c, len(pts), (Q0, Q1, Qh + Qh2, P(0), P(1))
for nu_s, R_s in [('3/25', '177/20'), ('1249/10000', '173/20')]:
    res = [comp(dps, N, nu_s, R_s) for dps, N in [(50, 400), (70, 800)]]
    (c1, n1, chk), (c2, n2, _) = res
    mp.dps = 70
    print(f"nu={nu_s} R={R_s}: c = {mp.nstr(c2, 30)} (nodes {n1} vs {n2}: diff {mp.nstr(abs(c1 - c2), 3)}); "
          f"checks Q(0)={mp.nstr(chk[0],12)} Q(1)={mp.nstr(chk[1],5)} Q(1/3)+Q(2/3)={mp.nstr(chk[2],12)} P(0)={chk[3]} P(1)={mp.nstr(chk[4],12)}")
    from flint import arb, ctx
    ctx.prec = 600
    rad = max(abs(c1 - c2) * 10, mpf(10) ** -45)
    cA = arb(mp.nstr(c2, 60)) + arb(0, mp.nstr(rad, 5))
    Rr = arb(int(R_s.split('/')[0])) / int(R_s.split('/')[1])
    kap = 1 - cA.log() / Rr
    th = arb(1) / 2; s2 = arb(2).sqrt(); x = th / s2
    cth = 2 - th / 2 - (x.cos() / x.sin()) / s2
    ks = 1 - 2 / (2 - cth)
    h = (3 + kap - ((1 - kap) * (9 - kap - 8 * cth)).sqrt()) / 4
    dh = ((1 - kap) / (9 - kap - 8 * cth)).sqrt() * (x.cos() / x.sin()) ** 2 / 2
    print(f"  Arb: kappa = {kap.str(25)}\n  c(1/2) = {cth.str(25)}, k*(1/2) = {ks.str(25)}\n  kappa - k* = {(kap - ks).str(15)}\n  h(1/2;kappa) = {h.str(22)}, dh/dtheta(1/2) = {dh.str(15)}")
    for t in ['501/1000', '51/100', '527/1000', '55/100', '6/10', '75/100']:
        tt = arb(int(t.split('/')[0])) / int(t.split('/')[1]); xx = tt / s2
        ct = 2 - tt / 2 - (xx.cos() / xx.sin()) / s2
        ht = (3 + kap - ((1 - kap) * (9 - kap - 8 * ct)).sqrt()) / 4
        print(f"    theta={t}: c={ct.str(10)} h={ht.str(10)}")
