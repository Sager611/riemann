# Independent 50-digit implementation of Bellotti's constants from the TeX formulas (arXiv 2405.12545v2), not from his code.
# usage: python3 indep.py tables.json <rows: 1,33 | all> out.json [modes, default code,flaws2,all]
# Fix names (a mode is a set of them):
#   c3two  (4.11): representatives on both sides of gamma_r: c3 -> 2 c3 (valid since H_floor(1/nu) <= 1.12 log(1/nu), checked)
#   Grig   l.933: |G(1/2+it)| <= sum|theta_j theta_k|[j,k]^(-1/2) <= Gc W/log^2 W, Gc = c^2 P <= 9.7572 (sqfree.py, prodP.py)
#   vw4    (3.10): residue term carries log^4(VW) = (v+w)^4 log^4 T, printed and coded (u+w)^4
#   d3L    (specialcase): d3 = 4(1+1.0061(u L+1.333+...))/((v-u)^2 L); the code divides by L^2 more
#   b1r2   Corollary l.211 with r = sqrt2 delta: the r^3 terms of b1 carry sqrt2 (dropped in the code)
#   d1VW   fracd2 must evaluate d1 at VW (code: contrd2(u,v,w,x) evaluates it at T^(w+x))
#   e6g    (3.10): e^{1/(6 gamma)} <= e^{1/(6*0.447 L)} (code: e^{1/(6L)})
#   d2ps   Lemma d2overn: partial summation with the u-dependent Cully-Hugill-Trudgian bound (the paper uses d1(X) for all u <= X)
#   c1tex  c1: TeX's 44.66 (=280.594/2pi) instead of the code's I2I3/(2pi) (both are the paper's; code value used unless set)
import sys, json
from mpmath import mp, mpf, sqrt, log, exp, pi, quad, inf, harmonic, floor
mp.dps = 50
MODES = {'code': set(), 'flaws2': {'c3two', 'Grig'},
         'all': {'c3two', 'Grig', 'vw4', 'd3L', 'b1r2', 'd1VW', 'e6g', 'd2ps'}}
def consts(T0, T1, a0, u, x, v, w, F, Gc=mpf('9.7572')):
    L = log(T0); s = mpf(a0); u, x, v, w = map(mpf, (u, x, v, w))
    def nu(T):
        lT = log(T)
        if T0 < exp(mpf('46.2')): return 1 / (mpf('5.558691') * lT)
        if T0 < exp(mpf('170.2')):
            J = lT / 6 + log(lT) + log(mpf('0.618'))
            return (mpf('0.04962') - mpf('0.0196') / (J + mpf('1.15'))) / (J + mpf('0.685') + mpf('0.155') * log(lT))
        if T0 < exp(481958): return log(lT) / (mpf('21.233') * lT)
        return 1 / (mpf('53.989') * lT**(mpf(2)/3) * log(lT)**(mpf(1)/3))
    r = sqrt(2) * (1 - s); k = sqrt(2) if 'b1r2' in F else 1
    b2 = 2*r*mpf('4.7908')/4 + 2*r*log(1 + 1/r)*(4/pi - 1)/(1 + 2*r) + 2*r*mpf('2.6908') + 2 + 64*r**2/(1 + 2*r)**2
    b1 = (sqrt(2)/2*log(T0 + r) + k*32*r**2/(1 + 2*r)**2*log(T0 + r) - k*16*r**2/(1 + 2*r)**2) / L
    D = [1/pi**2, mpf('0.746'), mpf('0.825'), mpf('0.462')]
    def d1(y):
        ly = log(y); return (D[0]*ly**3 + D[1]*ly**2 + D[2]*ly + D[3] + mpf('9.73')*y**(-mpf(1)/4)*ly + mpf('0.73')*y**(-mpf(1)/2)) / ly**3
    VW = T0**(v + w); l = log(VW)
    d1VW = d1(VW) if 'd1VW' in F else d1(T0**(w + x))
    if 'd2ps' in F:   # A(X)/X + int_1^X A(u)/u^2 du, A(u) <= u(D1 l^3 + ...) + 9.73 u^(3/4) log u + 0.73 u^(1/2) (u >= 2), A = 1 on [1,2)
        P4 = lambda z: D[0]*z**4/4 + D[1]*z**3/3 + D[2]*z**2/2 + D[3]*z
        tail = mpf('9.73')*(4*2**(-mpf(1)/4)*log(2) + 16*2**(-mpf(1)/4)) + mpf('0.73')*2/sqrt(2)
        d2 = (d1(VW)*l**3 + mpf(1)/2 + P4(l) - P4(log(2)) + tail) / l**4
    else:
        d2 = d1VW / l + d1VW / 4
    d3 = 4*(1 + mpf('1.0061')*(u*L + mpf('1.333') + 11/sqrt(u*L*T0**u))) / ((v - u)**2 * L)
    if 'd3L' not in F: d3 = d3 / L**2
    q = v / u
    d4 = mpf('3.09')*(u + v)/(v - u)*(mpf('1.301')*(q**2 + 1) + mpf('1.084')*(q + 1) - mpf('0.116'))/(q - 1) + d3/(exp(T0**(u + v - x))*L) + d3/L*(2 + (x - (u + v))*L)
    d5 = mpf('1.0061')*(w*L + mpf('1.333') + mpf('3.95')/T0**(w/2)) / (w*L)
    G1 = mpf('0.632121') + quad(lambda y: y**(2 - 2*s)*exp(-y), [1, inf])
    a = 2*x*L*nu(T0); FA = exp(a - 1)/a if a < 1 else mpf(1)   # (4.9) monotonicity repair (NOTE 1), applied with c3two
    c2 = G1*d4*d5/w * ((1 + 5*(FA - 1)) if 'c3two' in F else 1)
    lnu1 = log(1/nu(T1))
    c3 = 2*G1*mpf('1.12')*d4*d5*lnu1/w * (2 if 'c3two' in F else 1)
    I = quad(lambda y: exp(1/(6*y))*y**(mpf('1.5') - 2*s)*exp(-pi*y/2), [1, inf])
    c4 = 2*sqrt(2*pi)*mpf('2.4')*d4*d5*I/w
    e6 = exp(1/(6*(2*s - mpf('1.5'))))
    J1 = sqrt(2*pi)*e6*mpf('66.7')*quad(lambda y: (2*y)**(mpf(27)/164)*y**(1 - 2*s)*exp(-pi*y/2), [3, inf]); J4 = J1
    J2 = mpf('1.461')*mpf('4.05206')
    J31 = sqrt(2*pi)*e6*mpf('66.7')*quad(lambda y: (y/3 + 1)**(mpf(27)/164)*y**(1 - 2*s)*exp(-pi*y/2), [1, inf])
    J32 = mpf('66.7')*mpf('1.0486')*mpf('3.62642')
    J3 = (2*J31 + J32)*T0**(mpf(27)/164)
    Gfac = Gc/(w*L) if 'Grig' in F else d5
    c5 = d4*Gfac/(2*pi)*T0**w/(w*L)*(T0**(x*(mpf('1.5') - 2*s)) + (T0**u/L**2)**(mpf('1.5') - 2*s))*(J1 + J2 + J3 + J4)
    fac = (v + w) if 'vw4' in F else (u + w)
    e6g = exp(1/(6*mpf('0.447')*L)) if 'e6g' in F else exp(1/(6*L))
    R1 = d2*fac**4*sqrt(2*pi)*e6g*L**4*T0**(x*(1 - s)) / (exp(mpf('0.447')*pi*L/2)*(mpf('0.447')*L)**(s - mpf('0.5')))
    if 'c1tex' in F: K = mpf('44.66')
    else:
        K = 2*(mpf('122.092') + exp(mpf(1)/6)*sqrt(2*pi)*mpf('66.7')*quad(lambda y: (1 + y/(mpf('0.447')*L))**(mpf(27)/164)*exp(-pi*y/2)*y**(-s), [1, inf]))/(2*pi)
    R2 = K*d2*(v + w)**4*T0**((v + w)/2 + x*(mpf('0.5') - s) + mpf(27)/164)*L**4
    c1 = exp(-L**2/T0**u) - R1 - R2 - mpf(10)**-30
    C = (c2 + c3 + c4)*(b1 + b2/(nu(T0)*L))/(c1**2 - c5)
    return dict(a=a, FA=FA, b1=b1, b2=b2, d2=d2, d3=d3, d4=d4, d5=d5, c1=c1, c2=c2, c3=c3, c4=c4, c5=c5, C=C, R1=R1, R2=R2,
                F4need=2 + 4*(1 - s)*exp(mpf(1)/6 - pi/2)/I, H=harmonic(floor(1/nu(T1))), Hb=mpf('1.12')*lnu1, cond4=c1**2 - c5)
tabs = json.load(open(sys.argv[1]))
rows = range(1, 39) if sys.argv[2] == 'all' else map(int, sys.argv[2].split(','))
modes = sys.argv[4].split(',') if len(sys.argv) > 4 else ['code', 'flaws2', 'all']
def T(sx):
    sx = sx.strip()
    return 3*mpf(10)**12 if sx.startswith('3') else exp(mpf(sx[5:-1]))
out = {}
for row in rows:
    cb, p = tabs['CB'][row - 1], tabs['uxvw'][row - 1]
    for m in modes:
        d = consts(T(cb[0]), T(cb[1]), cb[2], *p[2:6], F=MODES[m])
        out['%d_%s' % (row, m)] = {k: mp.nstr(vv, 25) for k, vv in d.items()}
        print(row, m, 'C =', mp.nstr(d['C'], 15), ' c1 %s c5 %s d4 %s R1 %s' % tuple(mp.nstr(d[k], 8) for k in ('c1', 'c5', 'd4', 'R1')), flush=True)
json.dump(out, open(sys.argv[3], 'w'), indent=1)
