#!/usr/bin/env python3
# certify2.py -- Level-2 certificate (PROOF.md section 5): ramp-smoothed cutoff of width lam (in log t) for the tail
# integral, plus order-k one-sided averages for the pointwise bound eta (Lemma 10_k). Library + endpoint search.
#   uv run --no-project --with python-flint --with mpmath python3 certify2.py [b0 ...]
import sys, time, math
from flint import arb, ctx
import certify as C
ctx.prec = 192
PI, LOG2, LOG4, LOG2PI = C.PI, C.LOG2, C.LOG4, C.LOG2PI
LN10 = arb(10).log()
LYB = C.LYB                       # log 1e19
LY1 = 38*LN10; LY2 = 76*LN10      # zones for the prime-power integral
H, TH, FO, FI = C.H, C.TH, C.FO, C.FI
T0 = C.T0

def Sm(T, m):         # >= sum_{|gamma|>T} |gamma|^-m  (m >= 2), from HSW
    L = T.log(); m1 = m - 1
    main = ((T/(2*PI)).log()/(m1*T**m1) + arb(1)/(m1**2*T**m1))/(2*PI)
    err = C.Rh(T)/T**m + C.CS*(m*L + 1)/(m*T**m)
    return 2*(main + err)

S1T0 = C.S1U(T0)
SM_T0 = {m: Sm(T0, m) for m in range(2, 8)}
_lowcache = {}
def lowz(k, dlt):     # S1(U) + (2/delta)^k S_{k+1}(U), U = 2/delta
    key = (k, dlt.str(20))
    if key not in _lowcache:
        U = 2/dlt
        _lowcache[key] = C.S1U(U) + (2/dlt)**k*Sm(U, k + 1)
    return _lowcache[key]

def eta_k(L, b0, k, dlt, e1=arb(0)):   # Lemma 10_k: |theta(x)-x| <= eta x ; 10 <= x <= 1e76 (e1 = eta_1 if x > 1e38)
    x = L.exp()
    zs = (L/2).exp()*(1 + k*dlt)**(H + k)*lowz(k, dlt) \
         + ((b0*L).exp() + ((1 - b0)*L).exp())/2*(1 + k*dlt)**(b0 + k)*(2/dlt)**k*SM_T0[k + 1]
    psi_err = k*dlt*x/2 + LOG2PI + 1/x**2 + zs
    psith = (1 + e1)*(x.sqrt() + x**TH + x**FO + L/LOG2*x**FI)
    return (psi_err + psith)/x

DGRID = [arb(10)**arb(e)/4 if False else (arb(10)**(arb(e)/4)) for e in range(-56, -19)]   # delta = 10^(e/4), 1e-14 .. 1e-5
def best_eta(L, b0, e1):
    best = None
    if L.upper() <= (LYB - arb('1e-3')).lower():  # Buthe 2018 (x <= 1e19), also on the window [x, x e^lam], lam <= 1e-3
        best = arb('1.95')/(L/2).exp()
    for k in (2, 3, 4, 5):
        for d in DGRID:
            if k*d > arb('1e-3'): continue
            v = eta_k(L, b0, k, d, e1)
            if best is None or v.upper() < best.upper(): best = v
    return best

def Mlow3(L, eta1):   # lower bound for M(x), 2.03e6 <= x <= 1e76 (M is decreasing in x)
    if not (L.upper() <= LY1.lower()) and not (L.lower() >= LY1.upper()):
        L = arb(L.upper())                       # box straddles 1e38: M(x) >= M(x_b)
    if L.upper() <= LY1.lower():
        return (C.Ia(H, L) - C.Ia(H, LY1)) - arb('1.95')*(C.Ia(FO, L) - C.Ia(FO, LY1)) + (1 - eta1)*(C.Ia(H, LY1) - C.Ia(H, LY2))
    if L.lower() >= LY1.upper():
        return (1 - eta1)*(C.Ia(H, L) - C.Ia(H, LY2))
    raise RuntimeError('Mlow3 box straddles 1e38')

def Mup3(L, eta1):    # upper bound for M(x), x <= 1e76 (M is decreasing in x)
    if not (L.upper() <= LY1.lower()) and not (L.lower() >= LY1.upper()):
        L = arb(L.lower())                       # box straddles 1e38: M(x) <= M(x_a)
    base = C.Ia(H, L) + C.Ia(TH, L) + C.Ia(FO, L) + C.Ja(FI, L)/LOG2
    far = (LOG4 - 1)*(C.Ia(H, LY2) + C.Ia(TH, LY2)) + LOG4*C.Ja(FO, LY2)/LOG2
    if L.upper() <= LY1.lower():
        mid = eta1*(C.Ia(H, LY1) + C.Ia(TH, LY1) + C.Ia(FO, LY1) + C.Ja(FI, LY1)/LOG2)
        return base + mid + far
    if L.lower() >= LY1.upper():
        return (1 + eta1)*base + far
    raise RuntimeError('Mup3 box straddles 1e38')

def nZphi(L, b0, lam):   # (sqrt x log x) * bound for |sum_rho F_phi(rho)/rho|
    kap = 1 + 2*(L + 2)/(L*(L + 1))
    low = (1 + 1/L)*(C.BETA0*(1 + kap/C.GAM1) + lam/2*S1T0)
    ratio = (L + 2)/(L*(L + 1))                    # |H'(L)|/H(L)
    high = (((b0 - H)*L).exp() + ((H - b0)*L).exp())/2*(1 + 1/L)*(2/lam + 3*ratio)*SM_T0[3]
    return low + high

def nWin(L, eta, lam):   # (sqrt x log x) * |window term| <= sqrt(x)(1+1/L) eta lam/2
    return (L/2).exp()*(1 + 1/L)*eta*lam/2

def lam_opt(L, b0, eta):  # float choice of the ramp width
    x12 = math.exp(float(L.mid())/2); a = math.exp((float(b0.mid()) - 0.5)*float(L.mid()))*float(SM_T0[3].mid())
    b = x12*float(eta.upper().mid())/2
    lam = math.sqrt(a/b) if b > 0 else 1e-6
    return min(max(lam, 1e-14), 1e-3)

def margins(kind, La, Lb, b0, eta1):
    L = C.box(La, Lb)
    LAMMAX = arb('1e-3')
    # e1 > 0 is needed as soon as the box or its window [x, x e^lam] (lam <= 1e-3) reaches beyond 1e38
    e1 = eta1 if (L.upper() + LAMMAX) > LY1.lower() else arb(0)
    eta = best_eta(L, b0, e1)
    if not (eta < arb('0.01')): return arb(-1)
    lam0 = lam_opt(L, b0, eta)
    best = None
    for f in (1.0, 0.5, 2.0):
        lam = arb(min(lam0*f, 1e-3))             # Lemma 11 is used with lam <= 1e-3 only
        if kind == 'N':
            nS2 = (L/2).exp()*L/(2*(L.exp() - 1)); ntiny = (-5*L/2).exp()
            m = (L/2).exp()*L*Mlow3(L + lam, eta1) - nZphi(L, b0, lam) - nWin(L, eta, lam) - nS2 - ntiny
        else:
            m = robin_pos(La, Lb, b0, eta, eta1) - (C.nMup_fn(L, eta1) if False else (L/2).exp()*L*Mup3(L, eta1)) \
                - nZphi(L, b0, lam) - LOG2PI*(-L/2).exp() - nWin(L, eta, lam) - C.nDup(L, eta)
        if best is None or m.lower() > best.lower(): best = m
    return best

def robin_pos(La, Lb, b0, eta, eta1):    # (sqrt x_a log x_a) * lower bound of sum_{x2<p<=x} p^-2 + L_2
    xa = La.exp(); xb = Lb.exp()
    x2lo, x2hi = C.x2_bounds(La, Lb)
    a = x2hi; YB = LYB.exp(); Y1 = LY1.exp()
    sp2 = arb(0)
    if a.upper() < YB.lower():                         # Buthe piece on (a, min(xa, 1e19)]
        b = xa if xa.upper() <= YB.lower() else YB
        lb = b.log()
        sp2 += C.E1(a.log()) - C.E1(lb) - arb('1.95')*b**(-arb(3)/2)/lb - arb('3.9')*a**(-arb(3)/2)/a.log()
        a2 = YB
    else:
        a2 = a
    if xa.lower() > YB.upper() and xa.lower() > a2.upper():   # eta_1 piece on (a2, min(xa, 1e38)]
        b = xa if xa.upper() <= Y1.lower() else Y1
        sp2 += C.E1(a2.log()) - C.E1(b.log()) - eta1*(1/(b*b.log()) + 4/(a2*a2.log()))
    th2lo = x2lo - arb('1.95')*x2lo.sqrt() if x2lo.upper() <= YB.lower() else x2lo*(1 - eta1)
    th1hi = xb if Lb.upper() <= LYB.lower() else xb*(1 + eta.upper())
    l2 = ((th2lo/th1hi).log1p()/th1hi.log()).log1p()
    return (La/2).exp()*La*(sp2 + l2)

def run2(kind, b0, L0, L1, eta1, w0=arb('0.25'), wmin=arb('1e-6')):
    La = L0; minm = None; nbox = 0; w = w0
    while La < L1:
        Lb = La + w
        if Lb > L1: Lb = L1
        m = margins(kind, La, Lb, b0, eta1); nbox += 1
        if m > 0:
            if minm is None or m.lower() < minm: minm = float(m.lower())
            La = Lb
            if w < w0: w = w*2
        else:
            w = w/2
            if w < wmin: return False, La, nbox, minm
    return True, L1, nbox, minm

def eta1_of(b0):     # sup of eta(t) over [1e19, 1e38] = eta(1e19) (decreasing), order 2, delta = 1e-9
    return eta_k(LYB, b0, 2, arb(10)**-9).upper()

if __name__ == '__main__':
    t0 = time.time()
    print('S3 =', SM_T0[3].str(6), 'S4 =', SM_T0[4].str(6), 'S5 =', SM_T0[5].str(6), 'S1(T0) =', S1T0.str(6))
    for b0s in (sys.argv[1:] or ['7/8']):
        b0 = arb(7)/8 if b0s == '7/8' else arb(b0s)
        e1 = eta1_of(b0)
        print(f'b0={b0s}: eta_1 = {e1.str(5)}', flush=True)
        for kind, L0 in (('N', 7*LN10), ('R', 10*LN10)):
            ok, Lend, nbox, minm = run2(kind, b0, L0, arb(170), e1)
            print(f'b0={b0s} {kind}: certified from 10^{float(L0/LN10):.0f} up to x = 10^{float(Lend/LN10):.6f} '
                  f'boxes={nbox} min={minm} t={time.time()-t0:.0f}s', flush=True)
