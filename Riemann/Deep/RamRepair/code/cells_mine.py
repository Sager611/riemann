# cells_mine.py -- own evaluation of the log-free zero-density bound (Padhi Thm 6.1 / Lemma 6.2, Bellotti's method) with the
# corrections of NOTE.md section 3:
#   N(1-delta_k, e^l) <= ntil * (d4 d5 / w) * kap2 * C / (c1^2 - c5) * X^{2 delta_k} / (delta_k l)
#   kap2 = (1 + (Y/X)^{2 nu})/2 (|X^w - Y^w| <= X^Re w + Y^Re w, Re w >= 2 nu), Y = U/l^2, nu = lam*(l)/l (own zfr table);
#   C = 1 (diagonal) + 2 delta Lp (near pairs, LP with repulsion d*) + 4 delta n(om) sum_i Gr(1 + i om) (far pairs, BOTH sides),
#       Gr(y) = |Gamma(1+2delta+iy)|/(Gamma(1+2delta) y) >= |Gamma(s+iy)| for 0 <= s <= 2 delta;
#   c5 uses |G(1/2+it)| <= 16 PG W/log^2 W (PG = prod_p (1+(sqrt p-1)/p^2)), c1 uses a corrected tail 2 e^{-l(l-3u)}/l^2.
import math, json, sys, os, glob, pickle
from flint import arb, acb, arb_series, ctx
ctx.prec = 100
PI = arb.pi(); GAM = arb.const_euler()
HERE = os.path.dirname(os.path.abspath(__file__))
# ---- zero-free region (own certified table) ----
_ZT = json.load(open(os.path.join(HERE, 'out', 'zfr_table.json')))
ZT = sorted((float(k), float(v)) for k, v in _ZT['table'].items()); ZABOVE = float(_ZT['above'])
def lam_at(l):
    for L, v in ZT:
        if L >= l - 1e-12: return v
    return ZABOVE
# ---- -zeta'/zeta(1+e) ----
_ZF = {}
def Zf(e):
    if e not in _ZF:
        s = arb_series([1 + arb(e), 1], prec=2).zeta(); _ZF[e] = -s[1]/s[0]
    return _ZF[e]
TINY = arb('1e-11')
def Mb(l, e):       # upper bound for sum_rho Re 1/(1+e+it-rho), |t| <= e^l (+ window slack): min of the two Hadamard bounds
    Z = Zf(e); lT = arb(l)
    m0 = (lT - (2*PI).log())/2 + Z
    Zm = float(Z.mid()); best = (float(m0.mid()), None)
    for j in range(50, 201):
        a = j/100; v = 0.5*(l - math.log(2*math.pi)) + (a*a + 0.5)/(2*a)*Zm + (l - math.log(math.pi))/(8*a)
        if v < best[0]: best = (v, a)
    if best[1] is not None:
        a = arb(best[1]); m1 = (lT - (2*PI).log())/2 + (a*a + arb('0.5'))/(2*a)*Z + (lT - PI.log())/(8*a)
        if m1.upper() < m0.upper(): m0 = m1
    return m0 + TINY
def cost_lo(d, e, de):
    a1 = arb(e); a2 = arb(e) + de; d = arb(d)
    v1 = a1/(a1*a1 + d*d); v2 = a2/(a2*a2 + d*d)
    return arb(min(float(v1.lower()), float(v2.lower())))
_GR = {}
def Gr(d, de):      # sup_{0<=s<=2de} |Gamma(1+s+id)| <= min(1, |Gamma(1+2de+id)|/Gamma(1+2de))  (upper bound, decreasing in |d|)
    if d <= 0: return arb(1)
    key = (d, de.mid().str(25))
    if key not in _GR:
        v = acb(1 + 2*de, arb(d)).gamma().abs_upper()/(1 + 2*de).gamma()
        v = arb(v.upper())
        _GR[key] = v if v < 1 else arb(1)
    return _GR[key]
def near_LP(h, de, l, e, dst):
    """upper bound for sum over representatives s != r with |gamma_s - gamma_r| <= 1 of |Gamma(w_rs)| X^{Re w_rs} / X^{2 de}
    (de = delta_k ball): LP over rectangles j >= 1 on both sides, weak duality."""
    B = Mb(l, e) - cost_lo(0, e, de)
    items = []
    for j in range(1, int(1/h) + 2):
        if (j + 1)*h <= dst: continue
        dm = max((j - 1)*h, dst)
        twod = float((2*de).upper())
        b = Gr(dm, de)/arb(dm) if dm >= twod else Gr(dm, de)/(2*de)
        items += [(b, cost_lo((j + 1)*h, e, de))]*2
    Bf = float(B.mid()); fl = sorted(((float(x.mid())/float(y.mid()), x, y) for x, y in items), key=lambda z: -z[0])
    mu = 0.0; Br = Bf
    for r, x, y in fl:
        if Br <= 0: break
        mu = r; Br -= float(y.mid())
    if Br > 0: mu = 0.0
    muA = arb(mu)
    ub = muA*(B if B > 0 else arb(0)) if mu > 0 else arb(0)
    for x, y in items:
        z = x - muA*y
        if float(z.upper()) > 0: ub += arb(z.upper())
    return ub
def ntil(h, de, l, e):
    return max(1, math.floor(float((Mb(l, e)/cost_lo(h/2, e, de)).upper())))
def nwin(om, de, l, e, h):
    d = arb(om)/2; a1 = arb(e); a2 = arb(e) + de
    f = min(float((a1/(a1*a1 + d*d)).lower()), float((a2/(a2*a2 + d*d)).lower()))
    n1 = math.floor(float((Mb(l, e)/arb(f)).upper()))
    n2 = math.ceil(om/h) + 1
    return min(n1, n2)
def far(om, de, l, e, h):
    """upper bound for sum over representatives s with |gamma_s - gamma_r| > 1 of |Gamma(w_rs)| X^{Re w_rs}/X^{2 de}:
    windows [1 + i om, 1 + (i+1) om] on BOTH sides, <= n(om) representatives each, Gr(y)/y decreasing."""
    n = nwin(om, de, l, e, h); acc = arb(0); i = 0
    while True:
        y = 1 + i*om
        if y >= 60: break
        acc += Gr(y, de)/arb(y); i += 1
    acc += (arb(-84)).exp()/(1 - (arb(-1.4)*om).exp())     # y >= 60: Gr(y)/y <= e^{-1.4 y}
    return 2*n*acc, n
_PG = None
def PG():           # prod_p (1 + (sqrt p - 1)/p^2), p <= 1e6 exactly, tail factor exp(2/sqrt(1e6))
    global _PG
    if _PG is None:
        N = 10**6; sv = bytearray([1])*(N + 1); sv[0] = sv[1] = 0
        for i in range(2, int(N**0.5) + 1):
            if sv[i]: sv[i*i::i] = bytearray(len(sv[i*i::i]))
        acc = arb(0)
        for p in range(2, N + 1):
            if sv[p]:
                pa = arb(p); acc += (1 + (pa.sqrt() - 1)/(pa*pa)).log()
        _PG = (acc + arb(2)/1000).exp()
    return _PG
IGAM = arb('5.922'); JGAM = arb('4.72')    # certified in gamma_cert.py (own)
def consts(p, de, l, d4fun):
    u, t, w, x = [arb(float(v)) for v in p]; v = t*u; al = 1 - de; L = arb(l)
    B1 = arb('10.45')*arb('0.618')*IGAM*L/(2*PI)*(L*((v + w)/2 + arb(1)/6 + x*(arb('0.5') - al))).exp()
    tail = 2*(-L*(L - 3*u)).exp()/(L*L) + arb('1e-100')
    c1 = (-(L*L)/(u*L).exp()).exp() - B1 - tail
    d4 = d4fun(u, t, x, L)
    d5 = arb('1.0061')*(1 + (arb('1.333') + 4*(-w*L/2).exp())/(w*L))
    Gh = 16*PG()*(w*L).exp()/(w*L)**2
    lY = u*L - 2*L.log()
    c5 = d4/(2*PI)*Gh*((L*x*(arb('1.5') - 2*al)).exp() + (lY*(arb('1.5') - 2*al)).exp())*arb('0.618')*(L/6).exp()*L*JGAM
    return dict(c1=c1, c5=c5, d4=d4, d5=d5, lY=lY, u=u, t=t, w=w, x=x)
def admissible(p, k, l, nu):
    u, t, w, x = [float(v) for v in p]; de = k/1e4
    return (de > nu and 0.02 < w < u and 1.3 <= t <= 2.5 and u*l >= 40 and u*l - 2*math.log(l) < x*l
            and 1.1*t*u <= x <= 2.6*t*u and 2*nu*x*l >= 1)
def dstar_none(l, k): return 0.0
def cell_parts(k, l, p, dstar=dstar_none, kp_list=(1.0, 1.5, 2.0, 2.5, 3.0, 4.0), om_list=(0.25, 0.35, 0.5)):
    """everything except d4: returns dict (arb strings) or None; logN = log(nt d4 d5 kap2 C/(w (c1^2 - d4 c5r))) + 2 x de l - log(de l)."""
    de_f = k/1e4; de = arb(k)/10**4; L = arb(l)
    nu_f = lam_at(l)/l
    if not admissible(p, k, l, nu_f): return None
    one = lambda u_, t_, x_, L_: arb(1)
    K = consts(p, de, l, one)          # d4 = 1: c5 below is c5/d4
    nu = arb(lam_at(l))/L
    kap2 = (1 + (((K['u'] - K['x'])*L - 2*L.log())*2*nu).exp())/2
    dst = dstar(l, k)
    best = None
    for kp in kp_list:
        h = kp*de_f
        en = min(((float(Mb(l, j*0.0005).mid())/float(cost_lo(h/2, j*0.0005, de).mid()), j*0.0005) for j in range(1, 400, 3)))[1]
        nt = ntil(h, de, l, en)
        el = min(((float(near_LP(h, de, l, j*0.004, dst).mid()), j*0.004) for j in range(1, 60, 2)))[1]
        Lp = near_LP(h, de, l, el, dst)
        fb = None
        for om in om_list:
            ef = min(((nwin(om, de, l, j*0.005, h), j*0.005) for j in range(1, 120, 2)))[1]
            fv, n = far(om, de, l, ef, h)
            if fb is None or fv.upper() < fb.upper(): fb = fv
        Cc = 1 + 2*de*Lp + 2*de*fb
        val = nt*Cc
        if best is None or val.upper() < best[0].upper(): best = (val, nt, Cc, kp)
    val, nt, Cc, kp = best
    st = lambda a: a.str(25, radius=True)
    return dict(k=k, l=l, nt=nt, C=st(Cc), kap2=st(kap2), d5=st(K['d5']), c1=st(K['c1']), c5r=st(K['c5']), x=repr(float(p[3])), w=repr(float(p[2])),
                u=repr(float(p[0])), t=repr(float(p[1])), kp=kp, dst=dst)
def logN(P, d4):
    """final upper bound for log N(1 - k/1e4, e^l) from cell_parts output P and an upper bound d4 (arb)."""
    de = arb(P['k'])/10**4; L = arb(P['l'])
    C, kap2, d5, c1, c5r = [arb(P[n]) for n in ('C', 'kap2', 'd5', 'c1', 'c5r')]
    x, w = arb(P['x']), arb(P['w'])
    den = c1*c1 - d4*c5r
    if not (c1 > 0 and den > 0): return None
    return (P['nt']*d4*d5*kap2*C/(w*den)).log() + 2*x*de*L - (de*L).log()
def cell(k, l, p, d4fun, dstar=dstar_none, kp_list=(1.0, 1.5, 2.0, 2.5, 3.0, 4.0), om_list=(0.25, 0.35, 0.5)):
    """upper bound (arb) for log N(1 - k/1e4, e^l); None if not admissible."""
    de_f = k/1e4; de = arb(k)/10**4; L = arb(l)
    nu_f = lam_at(l)/l
    if not admissible(p, k, l, nu_f): return None
    K = consts(p, de, l, d4fun)
    c1, c5 = K['c1'], K['c5']
    if not (c1 > 0 and c1*c1 - c5 > 0): return None
    nu = arb(lam_at(l))/L
    kap2 = (1 + (((K['u'] - K['x'])*L - 2*L.log())*2*nu).exp())/2      # (1 + (Y/X)^{2 nu})/2, log(Y/X) = (u - x) l - 2 log l
    dst = dstar(l, k)
    best = None
    for kp in kp_list:
        h = kp*de_f
        # float choices of the auxiliary epsilons (any choice is valid)
        en = min(((float(Mb(l, j*0.0005).mid())/float(cost_lo(h/2, j*0.0005, de).mid()), j*0.0005) for j in range(1, 400, 3)))[1]
        nt = ntil(h, de, l, en)
        el = min(((float(near_LP(h, de, l, j*0.004, dst).mid()), j*0.004) for j in range(1, 60, 2)))[1]
        Lp = near_LP(h, de, l, el, dst)
        fb = None
        for om in om_list:
            ef = min(((nwin(om, de, l, j*0.005, h), j*0.005) for j in range(1, 120, 2)))[1]
            fv, n = far(om, de, l, ef, h)
            if fb is None or fv.upper() < fb.upper(): fb = fv
        Cc = 1 + 2*de*Lp + 2*de*fb
        val = nt*Cc
        if best is None or val.upper() < best[0].upper(): best = (val, nt, Cc, kp)
    val, nt, Cc, kp = best
    lg = (nt*K['d4']*K['d5']*kap2*Cc/(K['w']*(c1*c1 - c5))).log() + 2*K['x']*de*L - (de*L).log()
    return lg, dict(nt=nt, C=Cc, kap2=kap2, d4=K['d4'], c5=c5, c1=c1, kp=kp)
