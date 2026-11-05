#!/usr/bin/env python3
# certify.py -- interval-arithmetic certificate for height-tied Nicolas and Robin ranges.
# Run: uv run --no-project --with python-flint --with mpmath python3 certify.py [b0 ...]
# Every inequality below is checked with arb balls (python-flint); see PROOF.md for the lemmas.
import sys, time
from flint import arb, ctx
import mpmath
ctx.prec = 192

def E1(z):            # exponential integral E_1(z) = -Ei(-z), z > 0
    return -((-z).ei())

PI = arb.pi(); EG = arb.const_euler(); E = arb.const_e()
LOG2 = arb(2).log(); LOG4 = arb(4).log(); LOG2PI = (2*PI).log()
T0 = arb(3) * arb(10)**12                 # Platt-Trudgian: RH for 0 < gamma <= 3e12
GAM1 = arb('14.1347')                      # lower bound for the least positive ordinate 14.134725...
BETA0 = 2 + EG - (4*PI).log()              # sum over all zeros of 1/(rho(1-rho)) = 0.0461914...
C1, C2, C3 = arb('0.1038'), arb('0.2573'), arb('9.3675')   # HSW 2022 N(T) bound, T >= e
CS = C1 + C2 + C3
LYB = arb(19) * arb(10).log()              # log 1e19 (Buthe 2018 range end)
LY = arb(38) * arb(10).log()               # log 1e38 (theta(sqrt t) range end)

def Rh(t):            # HSW: |N(t) - t/(2pi) log(t/(2 pi e))| <= Rh(t)
    L = t.log(); return C1*L + C2*L.log() + C3

def S2T(T):           # >= sum_{|gamma|>T} gamma^-2
    L = T.log()
    return ((T/(2*PI)).log() + 1)/(PI*T) + 2*Rh(T)/T**2 + CS*(2*L + 1)/T**2

def S3T(T):           # >= sum_{|gamma|>T} |gamma|^-3
    L = T.log()
    return (2*(T/(2*PI)).log() + 1)/(4*PI*T**2) + 2*Rh(T)/T**3 + 2*CS*(3*L + 1)/(3*T**3)

def S1U(U):           # >= sum_{0<|gamma|<=U} |gamma|^-1   (U >= 14)
    a = arb(14); c = 2*PI*E
    M = U/(2*PI)*(U/c).log()
    return 2*((M + Rh(U))/U + ((U/c).log()**2 - (a/c).log()**2)/(4*PI) + CS*(a.log() + 1)/a)

ST = S2T(T0); S3 = S3T(T0)

def Ia(a, L):         # int_x^inf t^(a-2) h(t) dt,  h = 1/log t + 1/log^2 t,  L = log x, a < 1
    c = 1 - a
    return a*E1(c*L) + (-c*L).exp()/L

def Ja(a, L):         # int_x^inf t^(a-2) h(t) log t dt
    c = 1 - a
    return (-c*L).exp()/c + E1(c*L)

H, TH, FO, FI = arb(1)/2, arb(1)/3, arb(1)/4, arb(1)/5

def nMlow(L):         # (sqrt x log x) * lower bound of M(x) = int_x^inf (psi-theta) w ; valid 2.03e6 <= x <= 1e38
    m = Ia(H, L) - Ia(H, LY) - arb('1.95')*(Ia(FO, L) - Ia(FO, LY))
    return (L/2).exp()*L*m

def nMup(L):          # (sqrt x log x) * upper bound of M(x) ; valid x <= 1e38
    m = (Ia(H, L) + Ia(TH, L) + Ia(FO, L) + Ja(FI, L)/LOG2
         + (LOG4 - 1)*(Ia(H, LY) + Ia(TH, LY)) + LOG4*Ja(FO, LY)/LOG2)
    return (L/2).exp()*L*m

def nZup(L, b0):      # (sqrt x log x) * bound for |sum_rho F(rho)/rho|
    kap = 1 + 2*(L + 2)/(L*(L + 1))
    low = (1 + 1/L)*BETA0*(1 + kap/GAM1)
    high = (((b0 - H)*L).exp() + ((H - b0)*L).exp())/2*(1 + 1/L)*(1 + kap/T0)*ST
    return low + high

# ---- Lemma C: pointwise |theta(x)-x| <= eta(x) x for x >= 1e19 (second-order averaging, delta = 1e-9)
DEL = arb(10)**-9
def setdelta(d):
    global DEL, LOWZ
    DEL = d; U = 2/DEL
    LOWZ = S1U(U) + 4/DEL**2*S3T(U)
setdelta(arb(10)**-9)
def eta_theta(L, b0):  # valid for 10 <= x <= 1e38 (used for x >= 1e19)
    x = L.exp()
    zs = (L/2).exp()*(1 + 2*DEL)**(arb(5)/2)*LOWZ \
         + ((b0*L).exp() + ((1 - b0)*L).exp())/2*(1 + 2*DEL)**(b0 + 2)*4/DEL**2*S3
    psi_err = DEL*x + LOG2PI + 1/x**2 + zs
    psith = x.sqrt() + x**TH + x**FO + L/LOG2*x**FI
    return (psi_err + psith)/x

def eta_of(L, b0):    # |theta(x) - x| <= eta x
    if (L.upper() <= LYB.lower()):
        return arb('1.95')/(L/2).exp()            # Buthe 2018, 1423 <= x <= 1e19
    return eta_theta(L, b0)

def nDup(L, eta):     # (sqrt x log x) * upper bound of D(x)
    return (L/2).exp()*eta**2*(1 + (1 + eta)**2/L)

# ---- CA numbers: F(y,k) = log(1 + 1/(y+...+y^k))/log y
def F1(y): return (1/y).log1p()/y.log()
def F2(y): return (1/(y + y*y)).log1p()/y.log()

def x2_approx(eps):   # float solution of F2(y) = eps
    mpmath.mp.dps = 40
    e = mpmath.mpf(eps)
    y = mpmath.sqrt(1/(e*20))
    for _ in range(60):
        y = mpmath.sqrt(1/(e*mpmath.log(y))) if y > 3 else y
        # refine with exact equation: log1p(1/(y+y^2)) = e log y
    f = lambda t: mpmath.log1p(1/(t + t*t)) - e*mpmath.log(t)
    y = mpmath.findroot(f, y)
    return y

def x2_bounds(La, Lb):  # rigorous x2lo <= x_2(F1(x_a)), x2hi >= x_2(F1(x_b))
    xa = La.exp(); xb = Lb.exp()
    ea = F1(xa); eb = F1(xb)
    out = []
    for (eps, side) in ((ea, 'lo'), (eb, 'hi')):
        y0 = x2_approx(float(eps.mid()) if False else mpmath.mpf(eps.mid().str(30, radius=False)))
        rel = mpmath.mpf('1e-15')
        for _ in range(40):
            yy = y0*(1 - rel) if side == 'lo' else y0*(1 + rel)
            ya = arb(mpmath.nstr(yy, 40))
            v = F2(ya)
            if side == 'lo' and v.lower() > eps.upper(): out.append(ya); break
            if side == 'hi' and v.upper() < eps.lower(): out.append(ya); break
            rel *= 4
        else:
            raise RuntimeError('x2 bound failed')
    return out[0], out[1]

# ---- margins on a box L in [La, Lb]
def box(La, Lb):
    return arb((La + Lb)/2, (Lb - La)/2)     # ball containing [La, Lb] (arb rounds outward)

def nicolas_margin(La, Lb, b0):
    L = box(La, Lb)
    nS2 = (L/2).exp()*L/(2*((L.exp()) - 1))           # S_2(x) <= 1/(2(x-1))
    ntiny = (L/2).exp()*(-3*L).exp()                   # |1/2 int log(1-t^-2) w| <= x^-3/log x
    return nMlow(L) - nZup(L, b0) - nS2 - ntiny

def robin_margin(La, Lb, b0):
    L = box(La, Lb)
    xa = La.exp(); xb = Lb.exp()
    x2lo, x2hi = x2_bounds(La, Lb)
    # sum_{x2 < p <= x} 1/p^2 lower bound (truncate at 1e19 if needed)
    if La.upper() <= LYB.lower():
        Lt = La
    elif La.lower() >= LYB.upper():
        Lt = LYB
    else:
        raise RuntimeError('box start overlaps log 1e19')
    sp2 = E1(x2hi.log()) - E1(Lt) - arb('1.95')*(-arb(3)/2*Lt).exp()/Lt \
          - arb('3.9')*x2lo**(-arb(3)/2)/x2lo.log()
    eta_box = eta_of(L, b0)                                          # |theta(x)-x| <= eta x on the box
    th1hi = xb*(1 + eta_box.upper()) if Lb.upper() > LYB.lower() else xb   # theta(x) < x for x <= 1e19
    th2lo = x2lo - arb('1.95')*x2lo.sqrt()
    l2 = ((th2lo/th1hi).log1p()/th1hi.log()).log1p()
    pos = (La/2).exp()*La*(sp2 + l2)
    neg = nMup(L) + nZup(L, b0) + LOG2PI*(-L/2).exp() + nDup(L, eta_box)
    return pos - neg

def run(kind, b0, L0, L1, w0=arb('0.25'), wmin=arb('1e-7')):
    f = nicolas_margin if kind == 'N' else robin_margin
    La = L0; minm = None; nbox = 0
    w = w0
    while La < L1:
        Lb = La + w
        if Lb > L1: Lb = L1
        m = f(La, Lb, b0); nbox += 1
        if m > 0:
            if minm is None or m.lower() < minm: minm = float(m.lower())
            La = Lb
            if w < w0: w = w*2
        else:
            w = w/2
            if w < wmin:
                return False, La, nbox, minm
    return True, L1, nbox, minm

def find_end(kind, b0, L0, Lguess):
    ok, Lfail, nbox, minm = run(kind, b0, L0, Lguess)
    return ok, Lfail, nbox, minm

if __name__ == '__main__':
    t0 = time.time()
    print('S_T =', ST.str(8), ' S3 =', S3.str(8), ' BETA0 =', BETA0.str(12), ' LOWZ(C) =', LOWZ.str(6))
    for b0s in (sys.argv[1:] or ['7/8']):
        setdelta(arb(10)**-9 if b0s != '1' else arb('1.25e-8'))
        if '/' in b0s:
            p, q = b0s.split('/'); b0 = arb(int(p))/int(q)
        else:
            b0 = arb(b0s)
        L0N = arb(7)*arb(10).log()            # Nicolas analytic range starts at 1e7
        L0R = arb(10)*arb(10).log()           # Robin analytic range starts at 1e10
        for kind, L0 in (('N', L0N), ('R', L0R)):
            ok, Lend, nbox, minm = run(kind, b0, L0, arb(80))
            print(f'b0={b0s} {kind}: certified from 10^{float(L0/arb(10).log()):.2f} up to log x = {Lend.str(10)} '
                  f'(x = 10^{float(Lend/arb(10).log()):.6f}), boxes={nbox}, min margin={minm}, t={time.time()-t0:.1f}s', flush=True)
