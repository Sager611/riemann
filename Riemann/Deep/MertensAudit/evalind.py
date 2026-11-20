# evalind.py -- independent re-evaluation (mpmath, 40 digits; no flint, no import of the audited code) of the MertensCH/HP141
# record constants from the DERIVATION.md formulas (Lemmas H1, H3-H5, CH piece) on the dumped float plan (parameters only).
import json, math, sys
from mpmath import mp, mpf, zeta, log, exp, pi, sqrt
mp.dps = 40
M78 = '/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78/'
P = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/MertensAudit/logs/plans.json'))
BP = mpf('0.87334585610671858'); D = mpf('0.875') - BP
TP = mpf(3000000000000); TPT = mpf(3000175332800); T1 = mpf(1000)
assert TP + mpf('0.52') < TPT
JL0 = mpf(P['JL0'].split('+/-')[0].strip('[ ')) + mpf('1e-3')        # their Arb upper end 19.5963, padded (checked separately: falsify.py)
# ---- RH strip (Lemma H1), sigma_0' = 1.02, eta' = 0.1, B = 1/2
kp = (mpf(1)/2 + mpf('0.1'))/2; Ap = log(3*zeta(mpf('1.1'))) + kp*log(mpf(2)/3); LZ = log(zeta(mpf('1.02'))/zeta(mpf('2.04')))
assert Ap > 0 and LZ > 0
sL = mpf('0.85')
def LamL(s): return (2*mpf('1.02') - s - mpf(1)/2)/(s - mpf(1)/2)
def logCL(s, sr): return (LamL(s) - 1)*Ap + LamL(s)*LZ + (sr + LamL(s)*((1 + kp)*(1 + mpf('1.02')) + 2))/T1
aL = (LamL(sL) - 1)*kp; CL = exp(logCL(sL, sL))
RHP = JL0 + 2*CL*TP**aL/aL
print('RH strip: Lam_L', mp.nstr(LamL(sL), 8), 'a_L', mp.nstr(aL, 8), 'log C_L', mp.nstr(log(CL), 8), 'log RHP', mp.nstr(log(RHP), 8))
# ---- half-plane certificates (Lemma A1 with B', Corollary A4 with T_1 := T'), translated by -D
def cert(fn, eta):
    d = json.load(open(M78 + fn)); eta = mpf(eta)
    s1 = mpf(d['sigma1']) - D
    ms = [(mpf(a) - D, mpf(b), mpf(c)) for a, b, c in d['masses']]
    assert min(m[0] for m in ms) > 1
    kap = (1 + eta - BP)/2; A = log(3*zeta(1 + eta)) - kap*log(1 + BP)
    Lam = sum(m[2] for m in ms); a = (Lam - 1)*kap
    R = s1; K = (Lam - 1)*A
    for sg, y, lam in ms:
        c = abs(y)
        if lam > 0: R += lam*((1 + kap)*(c + 1 + sg) + 2*c + 2); K += lam*log(zeta(sg)/zeta(2*sg))
        else: R += -lam*((3 + 2*kap)*c + sg); K += -lam*log(zeta(sg))
    assert TP >= 2*max(abs(m[1]) for m in ms)
    c2 = sum(m[2]*(m[0] - BP) for m in ms) - (s1 - BP); assert c2 > 0
    return dict(s=s1, a=a, Lam=Lam, lC=K + R/TP, lC1000=K + R/1000, kap=kap, A=A, CS=exp(logCL(sL, s1)))
C = {'ML093g5_e10_Bp': cert('c_ml093_g5.json', '0.1'), 'ML092g8_e05_Bp': cert('c_ml092_g8.json', '0.05')}
for k, v in C.items(): print(k, 's', mp.nstr(v['s'], 12), 'Lam', mp.nstr(v['Lam'], 10), 'a', mp.nstr(v['a'], 8), 'logCinf(T\')', mp.nstr(v['lC'], 10), 'logCinf(1000)', mp.nstr(v['lC1000'], 10))
beta = 3/(2*pi**2)
def X10(n): return mpf(10)**(mpf(n)/10**6)
def chdec(lch):
    l10 = lch/math.log(10); e = math.floor(l10); return mpf("%.6fe%d" % (10**(l10 - e), e))
def G(a, tau, T): return tau*T**(a - 1)/(1 - a) if tau <= T else (tau**a - T**a)/a + tau**a/(1 - a)
for th_s, p in P['plans'].items():
    th = mpf(th_s); names = P['names']
    n0 = int(round(p['L0']*10**6)); n1 = int(round(p['L1']*10**6)); X0 = X10(n0)
    cCH = max(3/(pi*mpf(10)**10) + mpf('11.39'), 3*X0**(1 - th)/(pi*mpf(10)**10) + mpf('11.39')*X0**(mpf(1)/2 - th))
    rg = p['ranges']; assert rg[0][0] == n0 and rg[-1][1] == n1 and all(rg[k+1][0] == rg[k][1] for k in range(len(rg)-1))
    worst = mpf(0); wk = None
    for k, (lo, hi, ci, lch) in enumerate(rg):
        ce = C[names[ci]]; s, a = ce['s'], ce['a']; ch = chdec(lch); xk, xk1 = X10(lo), X10(hi)
        assert lo < hi and ch <= xk**(1 - th) and xk >= 2*10**8 and sL <= th < 1
        r = 1 + (ch/2)*xk**(th - 1)
        dec = (xk if s <= th else xk1)**(s - th)*r**s
        tau1 = 2*xk1**(1 - th)/ch + 1
        ck = beta*ch + 6*xk**(mpf(1)/2 - th)*sqrt(r) + xk**(sL - th)*r**sL*RHP/(2*pi) \
             + (s - sL)*ce['CS']*TP**aL*dec*min(1/TP, tau1/TP**2)/pi + exp(ce['lC'])/pi*dec*G(a, tau1, TP)
        if ck > worst: worst, wk = ck, k
    nm, ci, lch, L1 = p['fin']; ce = C[nm]; s, a = ce['s'], ce['a']; ch = chdec(lch); X1 = X10(n1)
    assert s + a*(1 - th) - th <= 0 and s <= th and ch <= X1**(1 - th)
    r = 1 + (ch/2)*X1**(th - 1)
    cfin = beta*ch + 6*X1**(mpf(1)/2 - th)*sqrt(r) + X1**(sL - th)*r**sL*RHP/(2*pi) + (s - sL)*ce['CS']*TP**aL*X1**(s - th)*r**s/(pi*TP) \
           + exp(ce['lC'])/pi*(1/a + 1/(1 - a))*X1**(s + a*(1 - th) - th)*r**s*(2/ch + X1**(th - 1))**a
    c = max(cCH, worst, cfin)
    print('theta', th_s, 'CH', mp.nstr(cCH, 12), 'ranges max', mp.nstr(worst, 12), '(range', wk, 'of', len(rg), ')', 'final', mp.nstr(cfin, 12), '=> c', mp.nstr(c, 8))
