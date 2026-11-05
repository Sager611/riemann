#!/usr/bin/env python3
# certify3.py -- Level-2 and Level-3 certificates with the Chirre-Helfgott pointwise bound (RobinCH/REPORT.md).
#   uv run --no-project --with python-flint --with mpmath python3 certify3.py [level] [b0 ...]
# Changes against RobinHeight/certify2.py (copied here unchanged, used as a library):
#  * eta: |theta(x) - x| <= eta x is now min(Buthe 1.95/sqrt x [window below 1e19], eta_CH(x)) with
#    eta_CH(x) = pi/(3e12) + 113.67/sqrt x + psith(x)/x  (Chirre-Helfgott, arXiv 2512.15709, Cor. "cor:psiexpl",
#    plus our psi - theta bound); Lemmas 10 and 10_k are no longer used.
#  * eta_1 (theta on [1e19, 1e38], prime-power zone 1e38..1e76, Robin's p^-2 sum, theta(x_2)) = eta_CH(1e19).
#  * level 2: Lemma 11 (linear ramp) exactly as in certify2.py.
#  * level 3: Lemma 11_k (ramp = CDF of a sum of k uniforms on [0, lam/k]); zeros above T0 cost
#    x^(beta-1) [H (2k/lam)^k |gamma|^-(k+2) + |H'| (lam/2 |gamma|^-2 + 2 |gamma|^-3)].
#  * level 4: as level 3 with the ramp centred at x (support [x e^(-lam/2), x e^(lam/2)]), so the window costs
#    eta H(L') E|V - lam/2| = eta H(L') (lam/k) MAD_k (<= lam/sqrt(12k)) instead of eta H(L) lam/2, and the sharper D bound (Lemma 13).
import sys, time, math
from fractions import Fraction
from flint import arb, ctx
import certify as C
import certify2 as K
ctx.prec = 192
PI, LOG2, LOG2PI = C.PI, C.LOG2, C.LOG2PI
LN10 = K.LN10; LYB = C.LYB; LY1 = K.LY1
H, TH, FO, FI = C.H, C.TH, C.FO, C.FI
T0 = C.T0                                  # 3e12: RH below this height (Platt-Trudgian: 3,000,175,332,800)
TCH = arb(3)*arb(10)**12                   # Chirre-Helfgott: |psi(x) - x| <= pi x/TCH + CCH sqrt x for all x >= 1
CCH = arb('113.67')
LAMMAX = arb('1e-3')
KMAX = 60
KMAX4 = 120

_sm = {}
def Sm(m):                                 # >= sum_{|gamma| > T0} |gamma|^-m
    if m not in _sm: _sm[m] = K.Sm(T0, m)
    return _sm[m]
S1T0 = K.S1T0

def psith(L, e1):                          # psi(x) - theta(x) <= psith, 1 <= x <= 1e76 (e1 = eta_1 needed above 1e38)
    x = L.exp()
    return (1 + e1)*(x.sqrt() + x**TH + x**FO + L/LOG2*x**FI)

def eta_ch(L, e1):                         # |theta(x) - x| <= eta_ch x  (x >= 1e19 used); decreasing in x
    x = L.exp()
    return PI/TCH + CCH/(L/2).exp() + psith(L, e1)/x

def eta1_ch():                             # eta_1 >= |theta(y)/y - 1| on [1e19, 1e38]
    return eta_ch(LYB, arb(0)).upper()

def best_eta3(L, e1):                      # valid on the box L and on its window [x, x e^lam], lam <= 1e-3
    v = eta_ch(L, e1)
    if L.upper() <= (LYB - LAMMAX).lower():
        b = arb('1.95')/(L/2).exp()        # Buthe 2018, 1423 <= x <= 1e19
        if b.upper() < v.upper(): v = b
    return v

def nZlow(L, lam):                         # (sqrt x log x) * zeros with |gamma| <= T0, (B1)
    kap = 1 + 2*(L + 2)/(L*(L + 1))
    return (1 + 1/L)*(C.BETA0*(1 + kap/C.GAM1) + lam/2*S1T0)

def nZhigh(L, b0, lam, k):                 # (sqrt x log x) * zeros with |gamma| > T0; k = 0: Lemma 11 (B2); k >= 1: Lemma 11_k
    ratio = (L + 2)/(L*(L + 1))            # |H'(L)|/H(L)
    P = (((b0 - H)*L).exp() + ((H - b0)*L).exp())/2*(1 + 1/L)
    if k == 0:
        return P*(2/lam + 3*ratio)*Sm(3)
    return P*((2*k/lam)**k*Sm(k + 2) + ratio*(lam/2*Sm(2) + 2*Sm(3)))

def nWin(L, eta, lam):                     # (sqrt x log x) * |W| <= sqrt(x)(1 + 1/L) eta lam/2
    return (L/2).exp()*(1 + 1/L)*eta*lam/2

_logS = {}
def choose(L, b0, eta, level):             # float choice of (k, lam): minimise window + high zeros
    Lf = float(L.mid()); bf = float(b0.mid()); ef = float(eta.upper().mid())
    logA = Lf/2 + math.log(1 + 1/Lf) + math.log(ef/2)                     # window: A lam
    logP = math.log((math.exp((bf - 0.5)*Lf) + math.exp((0.5 - bf)*Lf))/2*(1 + 1/Lf))
    cands = []
    ks = [0] if level == 2 else range(1, KMAX + 1)
    for k in ks:
        if k == 0:
            logB = logP + math.log(2) + float(Sm(3).log().mid()); kk = 1
        else:
            if k not in _logS: _logS[k] = float(Sm(k + 2).log().mid())
            logB = logP + k*math.log(2*k) + _logS[k]; kk = k
        loglam = (math.log(kk) + logB - logA)/(kk + 1)
        loglam = min(loglam, math.log(1e-3))
        cost = math.exp(logA + loglam) + math.exp(logB - kk*loglam)
        cands.append((cost, k, math.exp(loglam)))
    cands.sort()
    return cands[:2]

def margins3(kind, La, Lb, b0, eta1, level):
    L = C.box(La, Lb)
    e1 = eta1 if (L.upper() + LAMMAX) > LY1.lower() else arb(0)
    eta = best_eta3(L, e1)
    if not (eta < arb('0.01')): return arb(-1), None
    best = None; arg = None
    for cost, k, lam0 in choose(L, b0, eta, level):
        for f in (1.0, 0.8, 1.25):
            lam = arb(min(lam0*f, 9e-4))    # lam < 1e-3 = LAMMAX strictly
            if kind == 'N':
                nS2 = (L/2).exp()*L/(2*(L.exp() - 1)); ntiny = (-5*L/2).exp()
                m = (L/2).exp()*L*K.Mlow3(L + lam, eta1) - nZlow(L, lam) - nZhigh(L, b0, lam, k) - nWin(L, eta, lam) - nS2 - ntiny
            else:
                m = K.robin_pos(La, Lb, b0, eta, eta1) - (L/2).exp()*L*K.Mup3(L, eta1) \
                    - nZlow(L, lam) - nZhigh(L, b0, lam, k) - LOG2PI*(-L/2).exp() - nWin(L, eta, lam) - C.nDup(L, eta)
            if best is None or m.lower() > best.lower(): best = m; arg = (k, float(lam.mid()))
    return best, arg

# ---------------- level 4: centred k-box ramp (Lemma 11_k with c = lam/2) + Lemma 13
_mad = {}
def MAD(k):     # E|S - k/2| for S = sum of k uniforms on [0,1] (Irwin-Hall), exact: (2/(k+1)!) sum_{j<=k/2} (-1)^j C(k,j) (k/2-j)^(k+1)
    if k not in _mad:
        v = Fraction(0)
        for j in range(k//2 + 1):
            v += (-1)**j*math.comb(k, j)*(Fraction(k, 2) - j)**(k + 1)
        v = 2*v/math.factorial(k + 1)
        assert 0 < v and v*v <= Fraction(k, 12)          # MAD <= sd = sqrt(k/12)
        _mad[k] = arb(v.numerator)/arb(v.denominator)    # exact rational -> arb ball
    return _mad[k]
def mcen(lam, k):   # m_c = E|V - lam/2| = (lam/k) MAD(k)  (V = (lam/k) S)
    return lam*MAD(k)/k
def Hf(u): return 1/u + 1/u**2
def Hpf(u): return 1/u**2 + 2/u**3          # |H'(u)|

def nZlow4(L, lam, k):                      # zeros with |gamma| <= T0, (B1') with m_c <= lam/sqrt(12k), x' = x e^(-lam/2)
    kap = 1 + 2*(L + 2)/(L*(L + 1)); c = lam/2; Lp = L - c
    return (1 + 1/L)*C.BETA0*(1 + kap/C.GAM1) + L*(c/2).exp()*Hf(Lp)*mcen(lam, k)*S1T0

def nZhigh4(L, b0, lam, k):                 # zeros with |gamma| > T0, (B2_k) at L' = L - lam/2, pairing in x'
    c = lam/2; Lp = L - c
    P = (L/2).exp()*L*(((b0 - 1)*Lp).exp() + (-b0*Lp).exp())/2
    return P*(Hf(Lp)*(2*k/lam)**k*Sm(k + 2) + Hpf(Lp)*(lam/2*Sm(2) + 2*Sm(3)))

def nWin4(L, eta, lam, k):                  # (sqrt x log x) eta H(L') m_c,  m_c = (lam/k) MAD(k) <= lam/sqrt(12k)
    Lp = L - lam/2
    return (L/2).exp()*L*eta*Hf(Lp)*mcen(lam, k)

def nDup2(L, eta):                          # (sqrt x log x) * Lemma 13 bound for D(x)
    return (L/2).exp()*eta**2/2*(1 + 2*eta/(3*(1 - eta)) + (1 + eta)**2/L + 2*eta*(1 + eta)**3/(3*L**2*(1 - eta*(1 + eta)/L)))

def choose4(L, b0, eta):
    Lf = float(L.mid()); bf = float(b0.mid()); ef = float(eta.upper().mid())
    logA0 = Lf/2 + math.log(1 + 1/Lf) + math.log(ef)
    logP = math.log((math.exp((bf - 0.5)*Lf) + math.exp((0.5 - bf)*Lf))/2*(1 + 1/Lf))
    cands = []
    for k in range(1, KMAX4 + 1):
        if k not in _logS: _logS[k] = float(Sm(k + 2).log().mid())
        logA = logA0 + math.log(float(MAD(k).mid())/k)
        logB = logP + k*math.log(2*k) + _logS[k]
        loglam = min((math.log(k) + logB - logA)/(k + 1), math.log(1e-3))
        cost = math.exp(logA + loglam) + math.exp(logB - k*loglam)
        cands.append((cost, k, math.exp(loglam)))
    cands.sort()
    return cands[:2]

def margins4(kind, La, Lb, b0, eta1):
    L = C.box(La, Lb)
    e1 = eta1 if (L.upper() + LAMMAX) > LY1.lower() else arb(0)
    eta0 = best_eta3(L, e1)
    if not (eta0 < arb('0.01')): return arb(-1), None
    best = None; arg = None
    for cost, k, lam0 in choose4(L, b0, eta0):
        for f in (1.0, 0.8, 1.25):
            lam = arb(min(lam0*f, 9e-4)); c = lam/2    # lam < 1e-3 = LAMMAX strictly
            eta = best_eta3(L - c, e1)        # covers x in the box and the window [x e^-c, x e^c] (eta decreasing)
            if kind == 'N':
                nS2 = (L/2).exp()*L/(2*(L.exp() - 1)); ntiny = (-5*L/2).exp()*(3*c).exp()*L/(L - c)
                m = (L/2).exp()*L*K.Mlow3(L + c, eta1) - nZlow4(L, lam, k) - nZhigh4(L, b0, lam, k) - nWin4(L, eta, lam, k) - nS2 - ntiny
            else:
                m = K.robin_pos(La, Lb, b0, eta, eta1) - (L/2).exp()*L*K.Mup3(L - c, eta1) - nZlow4(L, lam, k) \
                    - nZhigh4(L, b0, lam, k) - LOG2PI*(-L/2).exp()*c.exp()*L/(L - c) - nWin4(L, eta, lam, k) - nDup2(L, eta)
            if best is None or m.lower() > best.lower(): best = m; arg = (k, float(lam.mid()))
    return best, arg

def run3(kind, b0, L0, L1, eta1, level, w0=arb('0.25'), wmin=arb('1e-6')):
    La = L0; minm = None; nbox = 0; w = w0; kmax_used = 0
    while La < L1:
        Lb = La + w
        if Lb > L1: Lb = L1
        m, arg = margins4(kind, La, Lb, b0, eta1) if level == 4 else margins3(kind, La, Lb, b0, eta1, level); nbox += 1
        if m > 0:
            if minm is None or m.lower() < minm: minm = float(m.lower())
            if arg: kmax_used = max(kmax_used, arg[0])
            La = Lb
            if w < w0: w = w*2
        else:
            w = w/2
            if w < wmin: return False, La, nbox, minm, kmax_used
    return True, L1, nbox, minm, kmax_used

def parse_b0(s):
    if '/' in s:
        p, q = s.split('/'); return arb(int(p))/int(q)
    return arb(s)

if __name__ == '__main__':
    t0 = time.time()
    level = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    e1 = eta1_ch()
    print('eta_1 (CH, 1e19) =', e1.str(6), ' eta_CH(1e40) =', eta_ch(40*LN10, e1).str(6), flush=True)
    for b0s in (sys.argv[2:] or ['0.87338196661320752', '7/8', '1']):
        b0 = parse_b0(b0s)
        for kind, L0 in (('N', 7*LN10), ('R', 10*LN10)):
            ok, Lend, nbox, minm, kmu = run3(kind, b0, L0, arb(170), e1, level)
            print(f'level {level} b0={b0s} {kind}: certified from 10^{float(L0/LN10):.0f} up to x = 10^{float(Lend/LN10):.6f} '
                  f'boxes={nbox} min={minm} kmax={kmu} t={time.time()-t0:.0f}s', flush=True)
