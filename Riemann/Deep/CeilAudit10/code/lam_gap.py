# lam_gap.py — ceilaudit10 T1: R(lam) = m1^2/(m0 m2) vs 2nu(lam) for zeta on dense grids in the gaps of explore63's C1.
# m_n(lam) = int_0^inf omega(e^{2u}) e^{u/2} cosh(u)^(-2lam) tanh(u)^(2n) du  (explore63 mp_lambda.py convention; Newton k=1 needs R >= 2nu N/(N-1)).
import sys, mpmath as mp
sys.path.insert(0, 'Riemann/Deep/ExplorePositivity')
from mp_core import omega
mp.mp.dps = 30
def moments(lam):
    s = 1/mp.sqrt(lam)
    pts = sorted(set([mp.mpf(0)] + [min(mp.mpf(3.3), s*k) for k in (0.25,0.5,1,1.5,2,3,4,6,8,12,16,24)] + [mp.mpf(v) for v in (0.25,0.5,0.75,1,1.5,2,2.5,3.3)]))
    pts = [p for p in pts if p <= 3.3]
    f = lambda u, n: omega(mp.exp(2*u)) * mp.exp(u/2) * mp.cosh(u)**(-2*lam) * mp.tanh(u)**(2*n)
    return [mp.quad(lambda u: f(u, n), pts) for n in range(3)]
def two_nu(lam):
    return (2*lam+2)*(2*lam+3)/(6*lam*(2*lam+1))
def row(lam):
    m = moments(lam); R = m[1]**2/(m[0]*m[2]); tn = two_nu(lam)
    return R, tn
mode = sys.argv[1]
if mode == 'check':
    for lam in (1.25, 30):
        R, tn = row(mp.mpf(lam)); print('check lam=%s R=%s 2nu=%s' % (lam, mp.nstr(R, 6), mp.nstr(tn, 6)), flush=True)
elif mode == 'low':
    worst = None
    for i in range(0, 26):
        lam = mp.mpf('1.10') + mp.mpf(i)/125
        R, tn = row(lam); mg = (tn - R)/tn
        if worst is None or mg < worst[0]: worst = (mg, lam, R, tn)
        if i % 5 == 0: print('lam=%s R=%s 2nu=%s margin=%s' % (mp.nstr(lam,5), mp.nstr(R,6), mp.nstr(tn,6), mp.nstr(mg,4)), flush=True)
    print('LOW [1.10,1.30] 26 pts: min relative margin (2nu-R)/2nu = %s at lam=%s (R=%s, 2nu=%s)' % (mp.nstr(worst[0],5), mp.nstr(worst[1],6), mp.nstr(worst[2],6), mp.nstr(worst[3],6)), flush=True)
else:
    worst = None; lam = mp.mpf(25); k = 0
    while lam <= 1e5:
        R, tn = row(lam); mg = (tn - R)/tn; sc = lam*(tn - R)
        if worst is None or mg < worst[0]: worst = (mg, lam)
        if k % 6 == 0: print('lam=%s R=%s 2nu=%s margin=%s lam*(2nu-R)=%s' % (mp.nstr(lam,6), mp.nstr(R,8), mp.nstr(tn,8), mp.nstr(mg,4), mp.nstr(sc,6)), flush=True)
        lam *= mp.mpf('1.15'); k += 1
    print('HIGH [25,1e5] %d pts (ratio 1.15): min relative margin = %s at lam=%s; last lam*(2nu-R) = %s (asymptotic 1/3)' % (k, mp.nstr(worst[0],5), mp.nstr(worst[1],6), mp.nstr(sc,6)), flush=True)
