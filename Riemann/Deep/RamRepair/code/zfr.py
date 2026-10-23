# zfr.py -- independent re-derivation (ball arithmetic, python-flint/Arb) of the per-height zero-free region
# (replaces Padhi's Thm 4.6) and of the two-zero repulsion margin (Thm 5.2), from the Bellotti-Trudgian-Yang (BTY) recipe.
# Every constant is recomputed here; see Riemann/Deep/RamRepair/NOTE.md section 3.1 for the derivation.
#   M1 = a1 (F(s-1+e) + F(s-e)) - [ e w0 ((a kap/2) l - G - c10) + a0 F(s-1) + Eterm + TAILS ],   e = eta = lam0/l, s = sigma
#   Eterm = Skap C1/(pi x0^2) (a lam0 e^2 + e^3 sum_{k>=1} a_k log(k/2)) + a0 Skap (C1 I0(x0)/(2 pi) + C2/sigma0^3) e^3
#   M1 > 0 on a box => no zero beta + i t with (1-beta) log t = lam0, log t = l in the box.
import math, json, sys, os, time
from fractions import Fraction as Fr
from flint import arb, acb, arb_series, ctx
from rig import up, dn, hi_max
ctx.prec = 80
PI = arb.pi()
th = arb(11338)/10000; kk = th.tan(); D = 2*th/kk; s2 = 1/th.cos()**2; thc = th/kk
def w_(u):    # u arb or arb_series
    return s2*(s2*(thc - u/2)*(u*kk).cos() + 2*thc - u + (2*th - u*kk).sin()/(2*th).sin() - 2 - 2*(th - u*kk).sin()/th.sin())
w0 = s2*(th*kk + 3*thc - 3)
KM = [Fr(1), Fr(-851, 859), Fr(780, 859), Fr(-525, 859), Fr(171, 859), Fr(28, 859), Fr(-29, 859)]
km = [arb(q.numerator)/q.denominator for q in KM]
kap = arb(433)/859; Skap = sum((abs(x) for x in km), arb(0))
CC = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
AKF = [Fr(sum(CC[j]*CC[j+k] for j in range(17-k))*(1 if k == 0 else 2), 14912370) for k in range(17)]
ak = [arb(q.numerator)/q.denominator for q in AKF]
a0, a1 = ak[0], ak[1]; a = sum(ak[1:], arb(0))
SLOGK2 = sum((ak[k]*(arb(k)/2).log() for k in range(1, 17)), arb(0))       # sum_{k>=1} a_k log(k/2)  (< 0)
K, T0 = 16, arb(10)**10
LOGH0 = arb(3*10**12).log()
TAILS = arb('1e-7')     # zeros |gamma| > K t + T0, pole terms k >= 1, far parts of E_k (k >= 1): certified by tails_bound() below
E2X = arb('0.035')      # extra eta^3 coefficient (per unit weight a_k) from U(Y) <= log(|Y|/2) + 0.35/|Y| (NOTE 3.1)

# ---------------- Laplace transform W(z) = int_0^D e^{-zu} w(u) du (closed form) ----------------
def J0(z): return -(-z*D).expm1()/z
def J1(z): return (1 - (-z*D).exp()*(1 + z*D))/(z*z)
IK = acb(0, kk)
def Sphi(z, phi): return (acb(0, phi).exp()*J0(z + IK) - acb(0, -phi).exp()*J0(z - IK))/acb(0, 2)
def W_ref(z):          # straightforward closed form (reference)
    z = acb(z)
    cosp = thc*(J0(z - IK) + J0(z + IK))/2 - (J1(z - IK) + J1(z + IK))/4
    return s2*(s2*cosp + (2*thc - 2)*J0(z) - J1(z) + Sphi(z, 2*th)/(2*th).sin() - 2*Sphi(z, th)/th.sin())
_E2M = acb(0, -2*th).exp(); _E2P = acb(0, 2*th).exp()       # e^{-ikD}, e^{+ikD}  (kD = 2 theta)
_EP1 = acb(0, th).exp(); _EM1 = acb(0, -th).exp(); _EP2 = _E2P; _EM2 = _E2M
_IS2 = 1/(2*th).sin(); _IS1 = 2/th.sin(); _I2 = acb(0, 2)
def W(z):              # same closed form, one exponential
    z = acb(z)
    E = (-z*D).exp(); Ep = E*_E2M; Em = E*_E2P
    zp = z + IK; zm = z - IK
    J0z = (1 - E)/z; J1z = (1 - E*(1 + z*D))/(z*z)
    J0p = (1 - Ep)/zp; J1p = (1 - Ep*(1 + zp*D))/(zp*zp)
    J0m = (1 - Em)/zm; J1m = (1 - Em*(1 + zm*D))/(zm*zm)
    cosp = thc*(J0m + J0p)/2 - (J1m + J1p)/4
    S2 = (_EP2*J0p - _EM2*J0m)/_I2; S1 = (_EP1*J0p - _EM1*J0m)/_I2
    return s2*(s2*cosp + (2*thc - 2)*J0z - J1z + S2*_IS2 - S1*_IS1)
def F(z, eta, sig):
    c = 2*sig - 1
    return sum((km[m]*Wf((z + c*m)/eta) for m in range(7)), acb(0))

# ---------------- constants ----------------
def _taylor(x, n=5):
    ctx.cap = n + 1
    return w_(arb_series([x, 1])).coeffs()
_c0 = _taylor(arb(0)); _cD = _taylor(D)
W2_0 = abs(2*_c0[2]); W3_0 = abs(6*_c0[3]); W3_D = abs(6*_cD[3])
W1_0 = _c0[1]
_I4 = {}
def I4(nu):            # upper bound of int_0^D |w''''(u)| e^{-nu u} du for Re z >= nu (nu double or ball; key = quarter grid <= nu)
    key = math.floor((dn(nu) if isinstance(nu, arb) else nu)*4)/4
    if key not in _I4:
        n = 400; Dm = float(D.upper())*(1 + 1e-15); acc = arb(0); nu_ = arb(key)
        for i in range(n):
            lo, hi = Dm*i/n, Dm*(i + 1)/n
            x = ball(lo, hi)
            sup = abs(24*_taylor(x, 4)[4]).upper()
            if key != 0: ex = ((-nu_*arb(lo)).exp() - (-nu_*arb(hi)).exp())/nu_
            else: ex = arb(hi - lo)
            acc += arb(sup)*arb(ex.upper())
        _I4[key] = arb(acc.upper())
    return _I4[key]
def Cw(nu, r):         # |W0(z)| <= Cw |z|^{-3}, Re z >= nu, |z| >= r  (four integrations by parts); nu, r doubles or balls (lower ends used)
    nl = dn(nu) if isinstance(nu, arb) else nu
    nu_ = arb(math.floor(nl*4)/4)
    rl = arb(r).lower()
    assert rl > 0
    return W2_0 + (W3_0 + (-nu_*D).exp()*W3_D + I4(nl))/rl
_I0 = {}
def _repsi(y):   # Re psi(1/4 + i y/2) at a point; increasing in |y| (series psi(a)+sum t^2/((n+a)((n+a)^2+t^2)))
    return acb(arb(1)/4, arb(y)/2).digamma().real
def qfloor(x, q=10**4):   # a double v on the grid Z/q with v <= every point of x (verified in Arb)
    k = math.floor((dn(x) if isinstance(x, arb) else x)*q) + 1
    while not (arb(k/q) <= arb(x).lower()): k -= 1
    return k/q
def I0(x0f):           # upper bound of int_R |Re psi(1/4 + iy/2)| (x0^2 + y^2)^{-3/2} dy, x0 >= x0f (decreasing in x0)
    key = qfloor(x0f)
    if key not in _I0:
        x0 = arb(key); acc = arb(0)
        def wint(lo, hi):   # exact int_lo^hi (x0^2+y^2)^{-3/2} dy
            F_ = lambda y: y/(x0*x0*(x0*x0 + y*y).sqrt())
            return F_(arb(hi)) - F_(arb(lo))
        grid = [0.001*i for i in range(3001)] + [3.0 + 0.01*i for i in range(1, 4701)] + [50.0 + i for i in range(1, 951)]
        vals = [abs(_repsi(y)) for y in grid]
        for i in range(len(grid) - 1):
            sup = vals[i] if vals[i] > vals[i+1] else vals[i+1]     # monotone Re psi => |Re psi| max at an endpoint
            acc += arb(sup.upper())*arb(wint(grid[i], grid[i+1]).upper())
        Y = arb(1000)
        tail = (2*Y.log() + 1)/(4*Y*Y)        # U(y) <= log(y/2) + 0.35/y <= log y for y >= 10 (NOTE 3.1);  int_Y^oo log(y) y^-3 dy
        _I0[key] = arb((2*(acc + tail)).upper())
    return _I0[key]
def G_const(sig0):     # sum_k a_k sum_m kappa_m(-log(pi)/2 + Re psi((s_k+delta_m)/2+1)/2) <= (a kap/2) log t - G
    s0 = arb(sig0)
    t1 = -a*kap/2*(2*PI).log() + kap/2*sum((ak[k]*arb(k).log() for k in range(1, 17)), arb(0))
    t0 = -kap*PI.log()/2
    for m in range(7):
        if KM[m] > 0: t0 += km[m]/2*(arb(m + 3)/2).digamma()
        else: t0 += km[m]/2*((s0 + (2*s0 - 1)*m)/2 + 1).digamma()
    tot = t1 + a0*t0 + arb('1e-10')        # + bound for the O(1/t) digamma errors (k >= 1)
    return -arb(tot.upper())               # G (rounded down)
# m_p: certified lower bounds of min_x G_p(1,x), G_p(1,x) = sum_{m<=15} p^-m P(mx); computed by mp_table() below
MP = None
def mp_table(path=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zfr_mp.json')):
    global MP
    if MP is None: MP = {int(k): arb(v).lower() for k, v in json.load(open(path)).items()}   # every point of the written ball is certified
    return MP
def phi_u(u, eta0, sig0):          # w(eta0 u) sum kappa_m e^{-d_m u} - kappa w0
    dm = [arb(m) if KM[m] > 0 else (2*arb(sig0) - 1)*m for m in range(7)]
    S = sum((km[m]*(-dm[m]*u).exp() for m in range(7)), arb(0))
    x = arb(eta0)*u
    if x > D: wv = arb(0)
    elif x < D: wv = w_(x)
    else: wv = arb(0, abs(w_(arb(D.mid(), x.rad() + D.rad()))).upper())
    return wv*S - kap*w0
def u0_cert(eta0, sig0, umax=80.0):
    """largest u0 (grid 0.05) with phi >= 0 on [0, u0], certified: near 0 via phi(0)=0, phi'(0) = -w0 sum kappa_m d_m > 0
    and a bound on |phi''|; elsewhere by bisection with ball evaluation."""
    import numpy as np
    dmf = [m if KM[m] > 0 else (2*sig0 - 1)*m for m in range(7)]
    us = np.linspace(0.0, umax, int(umax*20) + 1)
    def pf(u):
        x = eta0*u
        wv = float(w_(arb(x)).mid()) if x < float(D.mid()) else 0.0
        return wv*sum(float(KM[m])*math.exp(-dmf[m]*u) for m in range(7)) - float(kap.mid())*float(w0.mid())
    vals = [pf(u) for u in us]
    u0 = umax
    for i in range(1, len(us)):
        if vals[i] < 1e-12: u0 = float(us[i-1]) - 0.05; break
    # certify: [0, u1] by Taylor, [u1, u0] by bisection
    dm = [arb(m) if KM[m] > 0 else (2*arb(sig0) - 1)*m for m in range(7)]
    d1 = -w0*sum((km[m]*dm[m] for m in range(7)), arb(0))      # phi'(0) (w'(0) = 0)
    assert abs(W1_0) < arb('1e-15') and d1 > 0      # w'(0) = 0 exactly (w even autocorrelation)
    u1 = 0.05
    ctx.cap = 3
    U = arb(u1/2, u1/2)
    ser = None
    # phi(u) as a series at the ball U: second derivative enclosure
    us_ = arb_series([U, 1])
    x = arb(eta0)*us_
    S = sum((km[m]*(-dm[m]*us_).exp() for m in range(7)), arb_series([arb(0)]))
    ph = w_(x)*S
    g2 = 2*ph.coeffs()[2]
    if not (d1 - arb(u1)/2*abs(g2) > 0): return None
    stack = [(u1, u0)]
    while stack:
        lo, hi = stack.pop()
        if phi_u(ball(lo, hi), eta0, sig0) > 0: continue
        if hi - lo < 1e-6: return None
        mid = (lo + hi)/2; stack += [(lo, mid), (mid, hi)]
    return u0
_C10 = {}
def c10_const(eta0, sig0):
    key = (round(eta0, 6), round(sig0, 6))
    if key not in _C10:
        u0 = u0_cert(eta0, sig0)
        assert u0 is not None and u0 > 20
        MPT = mp_table(); s = arb(0)
        for p, m in MPT.items():
            if arb(15)*arb(p).log() <= arb(u0): s += arb(p).log()*m
        _C10[key] = (arb((kap*s).lower()), u0)
    return _C10[key][0]
def lemma5_ok(eta0, sig0):      # eps0 = Cw eta0^2/(w0 (2 sig0 - 1)^2) <= 1/500 (B(y) > 0 certified in lemma5 check)
    c = 2*arb(sig0) - 1; E0 = arb(eta0)
    e0 = Cw(c/E0, c/E0)*E0**2/(w0*c**2)
    return e0 < arb(1)/500
_CON = {}
def consts(eta0, sig0):
    """(E2a, E3, c10, G) valid for every eta <= eta0, sigma >= sig0  (eta0, sig0 already rounded outward)."""
    key = (eta0, sig0)
    if key not in _CON:
        assert eta0 <= 0.0145 and sig0 >= 0.993, (eta0, sig0)          # ranges of Lemma 5 and of tails_bound()
        S0 = arb(sig0); E0 = arb(eta0); x0 = S0 - arb(1)/2
        C1 = Cw(x0/E0, x0/E0); C2 = Cw(S0/E0, S0/E0)
        assert lemma5_ok(eta0, sig0), (eta0, sig0)
        E2 = Skap*C1/(PI*x0**2)                   # multiplies (a lam0 e^2 + (SLOGK2 + E2X a) e^3)
        E3 = a0*Skap*(C1*I0(x0)/(2*PI) + C2/S0**3)
        _CON[key] = (arb(E2.upper()), arb(E3.upper()), c10_const(eta0, sig0), G_const(sig0))
    return _CON[key]
def eta0_of(e):        # grid double (step 1e-4) >= every point of the ball e (verified in Arb)
    k = math.ceil(up(e)*10**4) - 1
    while not (arb(k/10**4) >= arb(e).upper()): k += 1
    return k/10**4
def sig0_of(sg):       # grid double (step 1e-4) <= every point of the ball sg
    return qfloor(sg)

W2s = 2*_c0[2]; W3s = 6*_c0[3]          # signed w''(0), w'''(0)
_P1 = {}
def P1(nu):     # upper bound of int_0^D u w(u) e^{-nu u} du  (>= |W'(z)| for Re z >= nu); w >= 0 decreasing on [0, D]
    key = math.floor(nu*4)/4          # nu: a double <= Re z (callers pass dn(...))
    if key not in _P1:
        n = 400; Dm = float(D.upper())*(1 + 1e-15); acc = arb(0); nu_ = arb(key)
        for i in range(n):
            lo, hi = Dm*i/n, Dm*(i + 1)/n
            # on [lo, hi]: u w(u) <= hi w(lo) if lo < D (w decreasing), and w = 0 beyond D; |w_(lo)| covers both cases
            sup = arb(hi)*arb(abs(w_(arb(lo))).upper()) if arb(lo) < D.upper() else arb(0)
            ex = ((-nu_*arb(lo)).exp() - (-nu_*arb(hi)).exp())/nu_ if key != 0 else arb(hi - lo)
            acc += arb(abs(sup).upper())*arb(ex.upper())
        _P1[key] = arb(acc.upper())
    return _P1[key]
def Wb(z):      # enclosure of W on a complex ball z: W(centre) + |W'| r  (mean value theorem)
    z = acb(z)
    rr = (z.real.rad()**2 + z.imag.rad()**2).sqrt()
    if rr == 0: return W(z)
    zc = acb(arb(z.real.mid()), arb(z.imag.mid()))
    nu = dn(z.real)
    e = arb(0, (arb(rr)*P1(nu)).upper())
    return W(zc) + acb(e, e)
_RW = {}
def Wf(z):
    """W(z); for Re z >= 100: W = w0/z + w''(0)/z^3 + w'''(0)/z^4 + O*((e^{-nu D}|w'''(D)| + I4(nu))/|z|^4), nu = floor(Re z)."""
    z = acb(z)
    nu = dn(z.real)
    if nu < 100: return Wb(z)
    nk = math.floor(nu)
    if nk not in _RW: _RW[nk] = arb(((-arb(nk)*D).exp()*W3_D + I4(nk)).upper())
    za = abs(z).lower()
    rr = _RW[nk]/arb(za)**4
    v = w0/z + W2s/z**3 + W3s/z**4
    rad = arb(0, rr.upper())
    return v + acb(rad, rad)

# ---------------- margins ----------------
def sigma_of(Lt, A):        # Lt arb (ball for log t), A arb
    T = K*Lt.exp() + T0
    return 1 - A/T.log()
def margin1(Lt, lam, A, C=None):
    e = lam/Lt; sig = sigma_of(Lt, A)
    if C is None:
        C = consts(eta0_of(e), sig0_of(sig))
    if not (A <= lam):                       # A <= lambda0 gives (1 - sigma)/eta < 1, used in tails_bound()
        raise ValueError('A > lambda0')
    E2, E3, c10, G = C
    lhs = e*w0*((a*kap/2)*Lt - G - c10) + a0*F(sig - 1, e, sig).real + E2*(a*lam*e*e + (SLOGK2 + E2X*a)*e**3) + E3*e**3 + TAILS
    rhs = a1*(F(sig - 1 + e, e, sig).real + F(sig - e, e, sig).real)
    return rhs - lhs
def margin2(Lt, lam0, lam1, d, A, C=None):
    """two zeros: primary beta0 + it (eta = lam0/l), second beta1 + i(t+d) with eta1 = lam1/l (same l-ball, widened)."""
    e = lam0/Lt; e1 = lam1/Lt; sig = sigma_of(Lt, A)
    if C is None:
        C = consts(eta0_of(e), sig0_of(sig))
    if not (A <= lam0): raise ValueError('A > lambda0')
    E2, E3, c10, G = C
    lhs = e*w0*((a*kap/2)*Lt - G - c10) + a0*F(sig - 1, e, sig).real + E2*(a*lam0*e*e + (SLOGK2 + E2X*a)*e**3) + E3*e**3 + TAILS
    rhs = a1*(F(sig - 1 + e, e, sig).real + F(sig - e, e, sig).real) \
        + a1*(F(acb(sig - 1 + e1, -d), e, sig).real + F(acb(sig - e1, -d), e, sig).real)
    return rhs - lhs
def ball(lo, hi):
    lo, hi = arb(lo), arb(hi); return (lo + hi)/2 + arb(0, ((hi - lo)/2).upper())

# ---------------- certified tail allowance (NOTE 3.1, review item 3) ----------------
def V2_const():        # int_0^D |w''(u)| du (upper sum, ball Taylor enclosures on 400 pieces)
    n = 400; Dm = float(D.upper())*(1 + 1e-15); acc = arb(0)
    for i in range(n):
        lo, hi = Dm*i/n, Dm*(i + 1)/n
        acc += arb(abs(2*_taylor(ball(lo, hi), 2)[2]).upper())*(arb(hi) - arb(lo))
    return arb(acc.upper())
def Umax():            # max over |Y| <= 10 of |Re psi(1/4 + iY/2)| (Re psi(1/4 + iy) increases with |y|)
    return arb(hi_max(abs(acb(arb(1)/4).digamma().real), abs(acb(arb(1)/4, arb(5)).digamma().real)))
def tails_bound(eta_cap=arb('0.0145'), lmax=46):
    """upper bound, uniform in eta <= eta_cap, H0 <= t <= e^lmax, sigma >= 0.993, A <= lambda0, for the terms covered by TAILS:
    (i) zeros with |gamma| > T = Kt + T0, all k; (ii) pole terms F(s_k - 1), k >= 1; (iii) parts |y| > kt - 10 of the E_k integrals and
    the terms W0(s/eta) of E_k, k >= 1.  Returns (total, parts)."""
    H0 = arb(3)*arb(10)**12; Tmax = K*arb(lmax).exp() + T0; e1 = arb(1).exp()
    V2 = V2_const(); eD = D.exp()
    # |Re W(x + iY)| <= Q/Y^2 for -1 < x < 7/eta (two integrations by parts of Re W = int e^{-xu} w(u) cos(Yu) du)
    Q = 4*w0*7/eta_cap + V2 + w0 + eD*(V2 + 2*w0 + D*w0)
    # zeros in (X, X + 1]: <= c2 log X + c1 (any explicit N(T) = (T/2pi) log(T/2pi e) + O*(E(T)), e.g. Hasanalizade-Shen-Wong 2022 or
    # Trudgian 2014 Cor. 1, gives 1/(2pi) log(X/2pi) + 2 E(X) <= log X + 40 for X >= 1e10)
    c1, c2 = arb(40), arb(1)
    def zsum(d0):      # >= sum_{gamma > T} (gamma - kt)^-2 with d0 = T - kt (unit intervals, log(T + 1 + j) <= log(T + 1) + log(1 + j/d0))
        Lg = (Tmax + 1).log()
        return (c1 + c2*Lg)*(1/d0 + 1/d0**2) + c2*(1/d0 + 1/(2*e1*d0**2))
    zer = arb(0)
    for k in range(K + 1):
        zer += ak[k]*2*zsum((K - k)*H0 + T0)      # factor 2: zeros below -T (distance |gamma| + kt is larger)
    zer *= Skap*Q*eta_cap**2
    pol = a*Skap*Q*eta_cap**2/H0**2                # |Re W((sigma - 1 + delta_m + ikt)/eta)| <= Q eta^2/(kt)^2
    C12 = Cw(arb(34), arb(34))                      # >= C1, C2 for all boxes (x0/eta0 >= 0.493/0.0145 > 34, Cw decreasing)
    R = H0 - 10
    efar = a*Skap*C12*eta_cap**3/(2*PI)*(Umax() + arb(3).log() + R.log() + arb(1)/2)/R**2
    t2 = a*Skap*C12*eta_cap**3/H0**3
    parts = {'zeros': zer, 'poles': pol, 'E_far': efar, 'W0_k>=1': t2}
    return zer + pol + efar + t2, parts
def G_tail():          # bound for the O(1/t) part of the archimedean terms, k >= 1 (Binet), covered by the 1e-10 in G_const
    H0 = arb(3)*arb(10)**12; X = arb(9)/2; Xl = arb('1.49'); Y = H0/2
    per = X/(2*Y**2) + 1/(24*Xl*Y) + X**2/(2*Y**2)  # |Re psi(z) - log(kt/2)| for 1.49 <= Re z <= 4.5, |Im z| >= H0/2
    return a*Skap*per/2

if __name__ == "__main__":
    # self-test: closed form W vs rigorous numerical integration
    for z in [acb('0.08'), acb('-0.9'), acb('135.2'), acb('0.3', '-7.1'), acb('140', '20')]:
        num = acb.integral(lambda u, an: (-z*u).exp()*w_(u), 0, D)
        print('W(%s): closed %s  integral %s' % (z.str(4), W(z).str(12), num.str(12)))
    print('w0', w0.str(15), ' w(0)', w_(arb(0)).str(15), " w'(0)", W1_0.str(5), " |w''(0)|", W2_0.str(12), " |w'''(0)|", W3_0.str(12), " |w'''(D)|", W3_D.str(12))
    print('a0 a1 a', a0.str(10), a1.str(10), a.str(10), ' kap', kap.str(10), ' Skap', Skap.str(10), ' sum a_k log(k/2)', SLOGK2.str(10))
    print('Cw(69.4,69.4)', Cw(69.4, 69.4).str(8), ' Cw(140,140)', Cw(140, 140).str(8))
    print('I0(0.4930)', I0(0.4930).str(8), ' I0(0.493516)', I0(0.493516).str(8))
    print('G(0.993)', G_const(0.993).str(8), ' G(0.9935164)', G_const(0.9935164).str(8))
    tb, parts = tails_bound()
    print('tails: total', tb.str(5), {k: v.str(3) for k, v in parts.items()}, ' allowance', TAILS.str(3), ' OK' if tb < TAILS else ' FAIL')
    gt = G_tail(); print('archimedean O(1/t) bound', gt.str(3), ' allowance 1e-10', ' OK' if gt < arb('1e-10') else ' FAIL')
    print('U_max', Umax().str(8), '  V2', V2_const().str(8))
    if not (tb < TAILS and gt < arb('1e-10')): raise SystemExit(1)
