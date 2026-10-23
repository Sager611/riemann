# d4_mine.py -- own upper bound for d4 = sum_n Psi(n)^2 n^{-1} (e^{-n/X} - e^{-n L^2/U})  (Padhi Thm 6.5, App. A; NOTE.md section 3.4):
#   d4 <= e^{-1/X} + D4 + REM + 2 sqrt(Sflat Ssharp) + Ssharp - e^{-L^2/U},  Sflat = e^{-1/X} + D4 + REM,
#   D4 = Delta^-2 [(Y2^2 - Y1^2)/2 - gamma Delta + k1 (Y1 + Y2 - 2 gamma) + 2 k0],  REM = R1 + ... + R8 (own constants, out/d4_consts.json),
#   k1, k0 upper ends from out/kconst_cert.json, M_I = 14.57, M_Z = 11.89 certified in b1_cert.py.
import math, json, os
from flint import arb, acb, ctx
ctx.prec = 100
HERE = os.path.dirname(os.path.abspath(__file__))
F = lambda v: arb(float(v))
eta, T2, W0, c, cp, s1, th1, b, th2, kap1 = [F(v) for v in (0.02, 400.0, 40.0, 0.02, 0.01, 0.01, 0.40, 0.38, 0.30, 0.38)]
PI = arb.pi(); GAM = arb.const_euler()
MI, MZ = arb('14.57'), arb('11.89')
RC = {k: arb(v) for k, v in json.load(open(os.path.join(HERE, 'out', 'd4_consts.json'))).items()}
KC = json.load(open(os.path.join(HERE, 'out', 'kconst_cert.json')))
K1UP, K0UP = arb(KC['k1_upper']), arb(KC['k0_upper'])
epsT = (-(eta*eta*T2*T2)/2).exp(); epsT2 = (-(eta*eta*(T2 - W0)**2)/2).exp()
sig2 = 1 + c - th1 - cp
Zb = MZ if MZ > 1/sig2 + 1/(1 - sig2) else 1/sig2 + 1/(1 - sig2)
gam0 = 4*arb('3.4')*W0.sqrt()*(-PI*W0/2).exp()/(PI*(1 - 1/(PI*W0)))     # int_{|v|>W0} |Gamma(sig+iv)| dv <= gam0*pi/2 ... (NOTE 4.5)
def Kfun(v):
    """|K_c(w)| for Re w = -th2, |Im w| = v, and for -th2 <= Re w <= c, |Im w| = W0 (v = W0)."""
    a_, a2 = kap1, kap1 - th2
    segc = (eta*eta*((kap1 + c)**2 + kap1**2)/2).exp()
    seg = PI*(a_ + a2)/(a_*a2*((a_ + a2)**2 + v*v))
    conn = 2*(kap1 + cp)*segc*epsT*epsT2/(T2**2*(T2 - W0)**2)
    return MZ/(2*PI)*(RC['AK']*MI*RC['iz1k_t2']*seg*segc + RC['AK2box']*MI*MI*conn)
_R = None
def Rconst(n=4000):
    """R = 2 int_0^W0 |Gamma(-th2 + iv)| (K(v) + e^{eta^2 th2^2/2}/th2^2) dv (upper Riemann sum: both factors decreasing in v)."""
    global _R
    if _R is None:
        Ew = (eta*eta*th2*th2/2).exp()/(th2*th2); h = W0/n; tot = arb(0)
        for i in range(n):
            v = arb(i)*h
            tot += arb(acb(-th2, v).gamma().abs_upper())*(Kfun(v) + Ew)
        _R = (arb((2*tot*h).upper()), arb((Kfun(W0) + (eta*eta*th2*th2/2).exp()/(W0*W0)).upper()))
    return _R
def remainder(lU, lV, lX):
    Dd = lV - lU; Y1, Y2 = lX - lV, lX - lU
    ex = lambda s, l: (s*l).exp()
    h = lambda sg: (ex(sg, lV) + ex(sg, lU))/Dd
    pi2, pi3 = (2*PI)**2, (2*PI)**3
    z1c = RC['z1c']
    out = {}
    out['R1'] = RC['gam_c_w1']*ex(c, lX)*z1c*Zb*MI*RC['iz1c_cp']*RC['A1']/pi3*h(-th1)*RC['J_mth1_0_T2_w1']*h(-cp)*RC['J_mcp_0_inf_w1']
    # R2: connectors Im s1 = +-T2 (pole of zeta(1+u) handled by |zeta(1+u)| <= 1/|u| + 1 + |1+u|/Re(1+u)) and |Im s1| > T2 on Re s1 = s1
    hbar = (ex(s1, lV) + ex(s1, lU))*(eta*eta*(th1*th1 - T2*T2)/2).exp()/(Dd*T2*T2)
    J0c = h(-cp)*RC['J_mcp_0_inf']; J1c = h(-cp)*RC['J_mcp_0_inf_w1']
    Zc = (th1 + s1)*((2 + (2 + W0 + T2)/sig2)*J0c + J1c/sig2) + 2*PI*((th1 + s1)**2 + 1).sqrt()*(ex(-cp, lV) + ex(-cp, lU))*(eta*eta*cp*cp/2).exp()/(Dd*cp*cp)
    out['R2'] = z1c*RC['iz1c_cp']*RC['A1box']/pi3*(2*RC['gam_c']*ex(c, lX)*MI*hbar*Zc
                + RC['gam_c_w1']*ex(c, lX)*h(s1)*RC['J_s1_T2_inf_w1']*RC['z1c_s1_cp']*RC['iz1c_s1']*J1c)
    out['R3'] = RC['gam_c']*ex(c, lX)*z1c*RC['Ax']*RC['iz1cp']*RC['iz1c_cp']/pi2*(ex(-(c - cp), lU) + ex(-(c - cp), lV))*(eta*eta*(c - cp)**2/2).exp()*epsT/(Dd*T2*T2)*J0c
    out['R4'] = 4*RC['gam_c']*ex(c, Y2)*z1c*RC['Ax']*RC['iz1cp']*RC['iz1c_cp']/(pi2*Dd*Dd)*2*(eta*eta*((c - cp)**2 + cp*cp)/2).exp()*epsT/(3*(T2 - W0)**3)
    rV = (-Dd*b).exp()*z1c*RC['AV']*MI*RC['iz1c_b']*(eta*eta*((c + b)**2 + b*b)/2).exp()/(2*b*(b + c)*(2*b + c))
    rVp = 2*(b + cp)*(cp*Dd).exp()*(eta*eta*((c + b)**2 + b*b)/2).exp()*epsT*epsT2/(2*PI*T2**2*(T2 - W0)**2)*z1c*RC['AVbox']*MI*RC['iz1c_cp']
    rU = (-Dd*b).exp()*z1c*RC['AU']*MI*RC['iz1b']*(eta*eta*((c - b)**2 + b*b)/2).exp()/(2*b*(b - c)*(2*b - c))
    rUp = 2*(b - cp)*(eta*eta*b*b).exp()*epsT*epsT2/(2*PI*T2**2*(T2 - W0)**2)*z1c*RC['AUbox']*RC['iz1cp']*MI
    out['R5'] = RC['gam_c']*ex(c, Y1)/(2*PI*Dd*Dd)*(rV + rVp)
    out['R6'] = RC['gam_c']*ex(c, Y2)/(2*PI*Dd*Dd)*(rU + rUp)
    Rr, Kw = Rconst()
    out['R7'] = (ex(-th2, Y1) + ex(-th2, Y2))*Rr/(2*PI*Dd*Dd)
    out['R8'] = gam0*ex(c, lX)*z1c**4/(2*PI*RC['z2_2c']) + 2*gam0*(th2 + c)*ex(c, Y2)*Kw/(Dd*Dd)
    return out
# ---------------- non-smooth part (Lemma A.3/A.4) ----------------
def _Phibar(x): return (x/arb(2).sqrt()).erfc()/2
def _phi(x): return (-(x*x)/2).exp()/(2*PI).sqrt()
def _varpi(x): return _phi(x) - x*_Phibar(x)                 # x >= 0
def _Phi1(y): return ((1 + y*y)*_Phibar(y) - y*_phi(y))/2    # int_y^oo varpi, y >= 0
C08 = arb('0.8')*(eta*eta/2).exp()
def Om_up(Alo, Ahi):
    """upper bound of Omega(A) = sum_r varpi((log r - A)/eta)/r for A in [Alo, Ahi] (Lemma A.3; bounds monotone in A)."""
    if Alo >= 0:
        return eta*(arb(1)/2 - _Phi1(arb(Ahi)/eta)) + C08*(-arb(Alo)).exp()
    if Ahi <= 0:
        y = arb(-Ahi)/eta
        return eta*_Phi1(y) + 2*_varpi(y)
    return eta/2 + C08
def Ssharp(lU, lV, lX, du=0.02):
    lUf, lVf, lXf = float(lU.mid()), float(lV.mid()), float(lX.mid())
    u0 = math.log(1e6)
    beta0 = eta + C08*((-(lU - arb(u0))).exp() + (-(lV - arb(u0))).exp())
    Sig = (arb(10**6).log() + arb('1.3326') + arb(4)/1000)*beta0**2
    u = u0; umax = lVf + 60*0.02
    while u < umax:
        ua, ub = u, u + du
        bet = arb(0)
        for K in (lU, lV):
            Kf = float(K.mid()); rk = float(K.rad()) + 1e-12
            bet += Om_up(Kf - ub - rk, Kf - ua + rk)
        Sig += bet**2*(arb(du) + 8*(-arb(ua)/2).exp())
        u = ub
    y0 = arb(umax - lVf)/eta
    Sig += 36*(eta + 2)**2*(-(y0*y0)).exp()/(2*PI)*(1 + eta/(2*y0))
    Dd = lV - lU
    best = None
    for k in range(1, 40):
        a = 0.25*k; e = a/lXf
        v = (e/math.e)**e*math.exp(a)*(1/e + 1)
        if best is None or v < best[0]: best = (v, a)
    a = arb(best[1]); e = a/lX
    return (e/arb(1).exp())**e*a.exp()*(1/e + 1)*(eta/Dd)**2*Sig + _Phibar(lU/eta)**2
def d4_bound(u, t, x, L, detail=False):
    """u, t (= tau), x arb; L arb (= l).  Returns arb upper bound of d4."""
    lU, lV, lX = u*L, t*u*L, x*L
    assert float(lU.lower()) >= 40 and lX > lV
    Y1, Y2, Dd = lX - lV, lX - lU, lV - lU
    D4m = ((Y2*Y2 - Y1*Y1)/2 - GAM*Dd + K1UP*(Y1 + Y2 - 2*GAM) + 2*K0UP)/(Dd*Dd)
    R = remainder(lU, lV, lX); rem = sum(R.values(), arb(0))
    Sf = (-1/lX.exp()).exp() + D4m + rem
    St = Ssharp(lU, lV, lX)
    d4 = Sf + 2*(Sf*St).sqrt() + St - (-(L*L)/lU.exp()).exp()
    if detail: return d4, dict(D4m=D4m, rem=rem, St=St, **R)
    return d4
