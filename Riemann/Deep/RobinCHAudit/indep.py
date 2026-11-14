#!/usr/bin/env python3
# indep.py (robinchaudit) -- independent re-implementation, from the formulas of RobinCH/PROOF.md and RobinHeight/PROOF.md
# (written from scratch, imports no RobinCH code), of RobinCH's Level-3 (centred k-box ramp, Lemma 11_k + Lemma 13) Nicolas
# and Robin margins on given boxes, in 256-bit arb (python-flint). Monotone terms at box endpoints, the rest on the ball.
#   uv run --no-project --with python-flint --with mpmath python3 indep.py boxes.json sel.json
import sys, json, math
from fractions import Fraction
from flint import arb, ctx
import mpmath
ctx.prec = 256

PI = arb.pi(); EG = arb.const_euler(); EE = arb(1).exp()
LN2 = arb(2).log(); LN4 = 2*LN2; LN10 = arb(10).log(); L2PI = (2*PI).log()
T0 = 3*arb(10)**12                      # H1: beta = 1/2 for |gamma| <= T0
G1 = arb('14.1347')                      # <= first ordinate 14.134725...
BETA0 = 2 + EG - (4*PI).log()            # H4: sum over zeros of 1/(rho(1-rho))
A1, A2, A3 = arb('0.1038'), arb('0.2573'), arb('9.3675'); AS = A1 + A2 + A3   # H3 (HSW), t >= e
LY19, LY38, LY76 = 19*LN10, 38*LN10, 76*LN10
Hh, Th, Qq, Fv = arb(1)/2, arb(1)/3, arb(1)/4, arb(1)/5
CHA, CHB = PI/(3*arb(10)**12), arb('113.67')   # H10: |psi(x) - x| <= CHA x + CHB sqrt x

def Rn(t): l = t.log(); return A1*l + A2*l.log() + A3
def Sm(m):      # >= sum_{|gamma|>T0} |gamma|^-m = 2 int_T0^inf t^-m dN; N = M + Q, |Q| <= Rn <= AS log t; m int R t^(-m-1) <= AS(log T + 1/m)/T^m
    T = T0; n = m - 1
    main = ((T/(2*PI)).log()/(n*T**n) + 1/(arb(n)**2*T**n))/(2*PI)
    return 2*(main + Rn(T)/T**m + AS*(T.log() + arb(1)/m)/T**m)
def S1():       # >= sum_{0<|gamma|<=T0} 1/|gamma| = 2[N(T0)/T0 + int_14^T0 N/t^2]
    U = T0; a = arb(14); c = 2*PI*EE
    return 2*((U/(2*PI)*(U/c).log() + Rn(U))/U + ((U/c).log()**2 - (a/c).log()**2)/(4*PI) + AS*(a.log() + 1)/a)
S1T0 = S1(); SM = {}
def sm(m):
    if m not in SM: SM[m] = Sm(m)
    return SM[m]

def E1(z): return -((-z).ei())
def Ia(a, L): c = 1 - a; return a*E1(c*L) + (-c*L).exp()/L          # int_x^inf t^(a-2) h dt
def Ja(a, L): c = 1 - a; return (-c*L).exp()/c + E1(c*L)             # int_x^inf t^(a-2) h log t dt
def P0(L): x = L.exp(); return x.sqrt() + x**Th + x**Qq + L/LN2*x**Fv   # psi - theta <= P0 (x <= 1e38)
def etach(L, fac): return CHA + CHB*(-L/2).exp() + fac*P0(L)/L.exp()  # Lemma 12
ETA1 = etach(LY19, arb(1)).upper()
def H(u): return 1/u + 1/u**2
def Hp(u): return 1/u**2 + 2/u**3

def MAD(k):     # E|S - k/2|, S Irwin-Hall(k): 2 int_0^{k/2} F_S, F_S(s) = (1/k!) sum_j (-1)^j C(k,j)(s-j)_+^k
    v = sum(Fraction((-1)**j*math.comb(k, j))*(Fraction(k, 2) - j)**(k + 1) for j in range(k//2 + 1))*2/math.factorial(k + 1)
    assert 0 < v and v*v <= Fraction(k, 12)
    return arb(v.numerator)/arb(v.denominator)

def Mlow(Ly):   # M(y) >= ..., point Ly, 2.03e6 <= y <= 1e76
    if Ly.upper() <= LY38.lower():
        return (Ia(Hh, Ly) - Ia(Hh, LY38)) - arb('1.95')*(Ia(Qq, Ly) - Ia(Qq, LY38)) + (1 - ETA1)*(Ia(Hh, LY38) - Ia(Hh, LY76))
    assert Ly.lower() >= LY38.upper()
    return (1 - ETA1)*(Ia(Hh, Ly) - Ia(Hh, LY76))
def Mup(Ly):    # M(y) <= ..., point Ly, y <= 1e76
    base = lambda L: Ia(Hh, L) + Ia(Th, L) + Ia(Qq, L) + Ja(Fv, L)/LN2
    far = (LN4 - 1)*(Ia(Hh, LY76) + Ia(Th, LY76)) + LN4*Ja(Qq, LY76)/LN2
    if Ly.upper() <= LY38.lower(): return base(Ly) + ETA1*base(LY38) + far
    assert Ly.lower() >= LY38.upper()
    return (1 + ETA1)*base(Ly) + far

def eta_box(La, Lb, c):  # >= |theta(t)/t - 1| on [x_a e^-c, x_b e^c]: eta_CH decreasing -> its value at La - c
    fac = 1 + ETA1 if (Lb + c).upper() > LY38.lower() else arb(1)
    v = etach(La - c, fac)
    if (Lb + c).upper() < LY19.lower():
        b = arb('1.95')*(-(La - c)/2).exp()     # Buthe, 1423 <= t <= 1e19
        if b.upper() < v.upper(): v = b
    return v

def F1(y): return (1/y).log1p()/y.log()
def F2(y): return (1/(y + y*y)).log1p()/y.log()
def x2_of(eps, side):    # side 'lo': y with F2(y) > eps (so y < x_2(eps)); 'hi': F2(y) < eps
    mpmath.mp.dps = 100
    e = mpmath.mpf(eps.mid().str(90, radius=False))
    y = mpmath.findroot(lambda t: mpmath.log1p(1/(t + t*t)) - e*mpmath.log(t), mpmath.sqrt(1/(e*mpmath.log(1/e)/2)))
    for rel in ('1e-60', '1e-50', '1e-40', '1e-30', '1e-20'):
        yy = y*(1 - mpmath.mpf(rel)) if side == 'lo' else y*(1 + mpmath.mpf(rel))
        ya = arb(mpmath.nstr(yy, 95))
        v = F2(ya)
        if side == 'lo' and v.lower() > eps.upper(): return ya
        if side == 'hi' and v.upper() < eps.lower(): return ya
    raise RuntimeError('x2 enclosure failed')

def box(La, Lb): return arb((La + Lb)/2, (Lb - La)/2)

def zeros_terms(L, b0, lam, k):          # normalised zero sums on the ball L: (B1') + Lemma 6 below T0, (B2_k) + pairing above
    c = lam/2; Lp = L - c; mc = lam*MAD(k)/k
    kap = 1 + 2*(L + 2)/(L*(L + 1))
    zlow = (1 + 1/L)*BETA0*(1 + kap/G1) + L*(c/2).exp()*H(Lp)*mc*S1T0
    pair = (L/2).exp()*L*(((b0 - 1)*Lp).exp() + (-b0*Lp).exp())/2
    zhigh = pair*(H(Lp)*(2*k/lam)**k*sm(k + 2) + Hp(Lp)*(lam/2*sm(2) + 2*sm(3)))
    return zlow, zhigh, mc

def nicolas(La, Lb, b0, k, lam):         # n(x)[M(x e^c) - |Z| - |W| - x'^-3/L' - 1/(2(x-1))] >= returned value on [La, Lb]
    c = lam/2; L = box(La, Lb); eta = eta_box(La, Lb, c)
    zlow, zhigh, mc = zeros_terms(L, b0, lam, k)
    M = (La/2).exp()*La*Mlow(Lb + c)    # n increasing, M decreasing
    win = (L/2).exp()*L*eta*H(L - c)*mc
    nS2 = (L/2).exp()*L/(2*(L.exp() - 1)); tiny = (L/2).exp()*L*(-3*(L - c)).exp()/(L - c)
    t = dict(M=M.lower(), zlow=zlow.upper(), zhigh=zhigh.upper(), win=win.upper(), eta=eta.upper())
    return M.lower() - (zlow + zhigh + win + nS2 + tiny).upper(), t

def robin(La, Lb, b0, k, lam):           # n(x_a)(sum p^-2 + L_2)_low - sup n(x)(M(x') + |Z| + log2pi/(x'L') + |W| + D)
    c = lam/2; L = box(La, Lb); eta = eta_box(La, Lb, c)
    xa, xb = La.exp(), Lb.exp()
    x2lo = x2_of(F1(xa), 'lo'); x2hi = x2_of(F1(xb), 'hi')
    YB, Y1 = LY19.exp(), LY38.exp(); a = x2hi; sp2 = arb(0)
    if a.upper() < YB.lower():           # Lemma 9 (Buthe) on (x2hi, min(xa, 1e19)]
        b = xa if xa.upper() <= YB.lower() else YB
        sp2 += E1(a.log()) - E1(b.log()) - arb('1.95')*b**(-arb(3)/2)/b.log() - arb('3.9')*a**(-arb(3)/2)/a.log()
        a2 = YB
    else: a2 = a
    if xa.lower() > YB.upper() and xa.lower() > a2.upper():   # eta_1 piece on (a2, min(xa, 1e38)]
        b = xa if xa.upper() <= Y1.lower() else Y1
        sp2 += E1(a2.log()) - E1(b.log()) - ETA1*(1/(b*b.log()) + 4/(a2*a2.log()))
    th2 = x2lo - arb('1.95')*x2lo.sqrt() if x2lo.upper() <= YB.lower() else x2lo*(1 - ETA1)
    th1 = xb if Lb.upper() <= LY19.lower() else xb*(1 + eta)
    l2 = ((th2/th1).log1p()/th1.log()).log1p()
    pos = (La/2).exp()*La*(sp2 + l2)
    zlow, zhigh, mc = zeros_terms(L, b0, lam, k)
    Mx = (Lb/2).exp()*Lb*Mup(La - c)
    lg = L2PI*(-L/2).exp()*c.exp()*L/(L - c)
    win = (L/2).exp()*L*eta*H(L - c)*mc
    D = (L/2).exp()*eta**2/2*(1 + 2*eta/(3*(1 - eta)) + (1 + eta)**2/L + 2*eta*(1 + eta)**3/(3*L**2*(1 - eta*(1 + eta)/L)))
    t = dict(pos=pos.lower(), M=Mx.upper(), zlow=zlow.upper(), zhigh=zhigh.upper(), win=win.upper(), D=D.upper(), eta=eta.upper())
    return pos.lower() - Mx.upper() - (zlow + zhigh + lg + win + D).upper(), t

def choose(kind, La, Lb, b0, kmax=120, fine=False):   # own float optimiser of window + high zeros, then rigorous margin
    Lf = float(((La + Lb)/2).mid()); bf = float(b0.mid())
    e = float(eta_box(La, Lb, arb(0)).upper().mid())
    logA0 = Lf/2 + math.log(Lf) + math.log(e) + math.log(1/Lf + 1/Lf**2)
    logP = Lf/2 + math.log(Lf) + math.log((math.exp((bf - 1)*Lf) + math.exp(-bf*Lf))/2) + math.log(1/Lf + 1/Lf**2)
    c = []
    for k in range(1, kmax + 1):
        logA = logA0 + math.log(float(MAD(k).mid())/k)
        logB = logP + k*math.log(2*k) + float(sm(k + 2).log().mid())
        ll = min((math.log(k) + logB - logA)/(k + 1), math.log(9e-4))
        c.append((math.exp(logA + ll) + math.exp(logB - k*ll), k, math.exp(ll)))
    c.sort(); best = None
    fs = [0.5 + 0.05*i for i in range(31)] if fine else [1.0, 0.8, 1.25]
    for _, k, l0 in c[:(4 if fine else 2)]:
        for f in fs:
            lam = arb(min(l0*f, 9e-4))
            m, t = (nicolas if kind == 'N' else robin)(La, Lb, b0, k, lam)
            if best is None or m > best[0]: best = (m, k, float(lam.mid()), t)
    return best


def A(me): return arb(me[0])*arb(2)**me[1]
def show(tag, kind, La, Lb, b0, k, lam, theirs=None):
    m, t = (nicolas if kind == 'N' else robin)(La, Lb, b0, k, arb(lam))
    lo = float((La/LN10).mid()); hi = float((Lb/LN10).mid())
    tt = ' '.join(f'{a}={float(v.mid()):.6g}' for a, v in t.items() if a in ('M', 'pos', 'zlow', 'zhigh', 'win', 'D', 'eta'))
    print(f'{tag} {kind} b0={float(b0.mid()):.5f} x in [10^{lo:.6f}, 10^{hi:.6f}] k={k} lam={lam:.4g}: ours m >= {float(m.lower()):.6g}'
          + (f'  (RobinCH m >= {theirs:.6g})' if theirs is not None else '') + f' {"PASS" if m > 0 else "FAIL"} | {tt}', flush=True)
    return m

if __name__ == '__main__':
    print(f'constants: S1(T0)={S1T0.str(8)} S2={sm(2).str(6)} S3={sm(3).str(6)} S42={sm(42).str(6)} BETA0={BETA0.str(12)} eta_1={ETA1.str(8)}')
    print(f'  E1 check: flint E1(25)={E1(arb(25)).str(15)} mpmath={mpmath.nstr(mpmath.e1(25), 15)}; E1(55.5)={E1(arb("55.5")).str(15)} mpmath={mpmath.nstr(mpmath.e1(mpmath.mpf("55.5")), 15)}')
    print(f'  MAD: k=1 {MAD(1)}, k=2 {MAD(2)}, k=40 {MAD(40).str(10)} vs normal approx {math.sqrt(40/12)*math.sqrt(2/math.pi):.10f}')
    D = json.load(open(sys.argv[1])); B = D['boxes']
    for s in D['summary']:
        rows = B[s['first']:s['last']]
        if s['X'] == 'probe':    # failing boxes at the endpoint: last 3 attempted
            for r in rows[-3:]:
                show('probe', r['kind'], A(r['La']), A(r['Lb']), arb(r['b0']), r['k'], float.fromhex(r['lam']), r['m'])
            continue
        ok = [r for r in rows if r['ok']]
        LN = math.log(10)
        pick = {'first': ok[0], 'last': ok[-1], 'min': min(ok, key=lambda r: r['m'])}
        for tgt in ((5e8 if s['kind'] == 'N' else None), 1e28):
            if tgt is None: continue
            for r in ok:
                if float(A(r['La']).mid()) <= math.log(tgt) <= float(A(r['Lb']).mid()): pick[f'x={tgt:.0e}'] = r
        print(f'== run {s["kind"]} b0={s["b0"]} to {s["X"]}: {len(ok)} passing boxes of {len(rows)} attempted')
        for tag, r in pick.items():
            show(tag, r['kind'], A(r['La']), A(r['Lb']), arb(r['b0']), r['k'], float.fromhex(r['lam']), r['m'])
    # Robin from 5e8 (as OddRobin's transfer uses RobinCH's margins there), own (k, lam) choice
    for b0 in (arb(1), arb('0.87338196661320752')):
        La = arb('20.03'); Lb = arb('20.04')
        m, k, lam, t = choose('R', La, Lb, b0)
        print(f'own-choice R b0={float(b0.mid()):.5f} log x in [20.03, 20.04]: k={k} lam={lam:.4g} m >= {float(m.lower()):.6g} (oddraudit F4: 0.468)')

