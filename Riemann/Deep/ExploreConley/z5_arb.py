# z5_arb.py — explore66: Arb (python-flint) certification of the Lehmer-pair block (z1/z2), rigorous ball arithmetic.
# (1) Taylor coefficients A_k of H_0(c+w)/scale from M rigorous samples on |w| = R, plus a rigorous aliasing bound from
#     M5 = max_{|w|=5} |H_0(c+w)| (cover of the circle by boxes); (2) truncation bound for k >= K under e^{-t d^2}, |t| <= smax, |w| <= r0;
# (3) rigorous: no zero of H_t, H_t' on |w| = r0 and empty exit set for all t in [-2e-3, 0]; zero counts 2 and 1 at t = 0;
# (4) collision window: H_{t1} takes a negative value on R and is positive at w = +-0.3 (real pair, t1 = t_c + 1e-9);
#     H_{t2} > 0 on [-0.3, 0.3] (t2 = t_c - 1e-9) => the pair is non-real at t2; hence t_c in (t2, t1).
import json, time, sys
from math import factorial
from flint import arb, acb, ctx
ctx.prec = 200
T0 = time.time()
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
zs = dict((n, s) for n, s in json.load(open(LOG + 'zeros_6700_6720.json')))
c = arb(zs[6709]) + arb(zs[6710])          # exact centre (a fixed ball of radius ~1e-60)
PI = acb(arb.pi()); I = acb(0, 1)
def H0(z):
    s = acb(0.5) + I*z/2
    return s*(s-1)*(PI**(-s/2))*(s/2).gamma()*s.zeta()/16
R = arb(2.5); M = 96; K = M//2
ws = [R*acb(arb(2*j)/M).exp_pi_i() for j in range(M)]
vals = [H0(c + w) for w in ws]
print('samples', round(time.time()-T0, 1), 's; max rel radius', max(float(v.rad()/abs(v).lower()) if abs(v).lower() > 0 else 9e9 for v in vals), flush=True)
A = []
for k in range(K):
    s = acb(0)
    for j in range(M): s += vals[j]*acb(arb(-2*j*k)/M).exp_pi_i()
    A.append(s/(M*R**k))
# rigorous M5: cover |w| = 5 by NB boxes
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
rho = arb(5); M5 = arb(0)
for j in range(NB):
    th = arb(2*j + 1)/NB                                  # centre angle (units of pi)
    cen = rho*acb(th).exp_pi_i()
    half = float(5*3.1416/NB)*1.05 + 1e-12               # >= half arc length (chord bound)
    box = acb(arb(cen.real.mid(), half), arb(cen.imag.mid(), half))
    v = abs(H0(c + box)).upper()
    if v > M5: M5 = arb(v)
print('M5 =', M5.str(5), ' time', round(time.time()-T0, 1), flush=True)
scale = A[2]
alias = M5*(R/rho)**M/(1 - (R/rho)**M)                    # |sum_{m>=1} a_{k+mM} R^{mM}| <= alias * rho^{-k} ... times R^0
A = [(A[k] + acb(arb(0, (alias*rho**(-k)).upper()), arb(0, (alias*rho**(-k)).upper())))/scale for k in range(K)]   # widen by aliasing, normalise
a = [x.real for x in A]                                    # H real on R: take real parts (imag parts are rigorous balls around 0)
print('max |Im A_k| upper:', max(float(abs(x.imag).upper()) for x in A), '  max rad a_k:', max(float(x.rad()) for x in a))
# truncation, rigorous for ALL k >= K. Coefficient bounds: A_geo(k) = (M10/|scale|) 10^{-k} (Cauchy, rigorous box-cover max M10
# on |w| = 10) and A_crude(k) = xi(1/2 + (c + rho')/2)/(8 |scale| rho'^k), rho' = k/1.75 (|xi(1/2+v)| <= xi(1/2+|v|): the Taylor
# coefficients of xi at 1/2 are >= 0). Heat majorant (Minkowski): |d^m e^{-t d^2} w^k| <= k!/(k-m)! (r0 + sig sqrt(k-m+1))^{k-m},
# sig = sqrt(2 smax), |w| <= r0, |t| <= smax.
smax = arb('0.0021'); r0 = arb('0.6'); sig = (2*smax).sqrt()
rho10 = arb(10); M10 = arb(0); NB10 = 6000
for j in range(NB10):
    cen = rho10*acb(arb(2*j + 1)/NB10).exp_pi_i(); half = float(10*3.1416/NB10)*1.05
    v = abs(H0(c + acb(arb(cen.real.mid(), half), arb(cen.imag.mid(), half)))).upper()
    if v > M10: M10 = arb(v)
sabs = abs(scale).lower()
print('M10 =', M10.str(5), ' M10/|scale| =', (M10/sabs).str(5), ' time', round(time.time()-T0, 1), flush=True)
def xi_real(x):     # rigorous xi(x) for real x > 1
    return x*(x-1)/2*(arb.pi()**(-x/2))*(x/2).gamma()*x.zeta()
def logterm_geo(k, m):
    return (M10/sabs).log() - k*rho10.log() + (arb(k+1).lgamma() - arb(k-m+1).lgamma()) + (k-m)*(r0 + sig*arb(k-m+1).sqrt()).log()
def logterm_crude(k, m):
    rp = arb(k)/arb('1.75')
    return xi_real(arb(1)/2 + (c + rp)/2).log() - arb(8).log() - sabs.log() - k*rp.log() + (arb(k+1).lgamma() - arb(k-m+1).lgamma()) + (k-m)*(r0 + sig*arb(k-m+1).sqrt()).log()
TR = []
for m in range(4):
    tot = arb(0); prev = None; mono = True
    for k in range(K, 3001):
        l = logterm_geo(k, m); tot += l.exp()
        if k > 200 and prev is not None and not (l < prev): mono = False
        prev = l
    g3001 = logterm_geo(3001, m)
    # k -> log geo_k is convex (second difference ~ sig/(2 sqrt(k)(r0 + sig sqrt k)) > 0), so on [3001, 16000] it is <= max of the endpoints
    g16000 = logterm_geo(16000, m)
    mid_sum = arb(13000)*(g3001.exp() if g3001 > g16000 else g16000.exp())
    c1_ = logterm_crude(16001, m); c2_ = logterm_crude(16002, m); q = (c2_ - c1_).exp()
    rat = [(logterm_crude(kk+1, m) - logterm_crude(kk, m)).exp().str(3) for kk in (16001, 32000, 64000, 128000)]
    far = c1_.exp()/(1 - q)
    T = tot + mid_sum + far
    TR.append(T)
    print('m=%d: tail <= %s  [sum_{K..3000} geo = %s (decreasing beyond k=200: %s); geo_3001 = e^%s, geo_16000 = e^%s; crude_16001 = e^%s, ratio %s, ratios at 16001/32000/64000/128000: %s]' % (
        m, T.str(4), tot.str(4), mono, g3001.str(5), g16000.str(5), c1_.str(5), q.str(3), rat), flush=True)
json.dump({'a_mid': [x.mid().str(40) for x in a], 'a_rad': [x.rad().str(5) for x in a], 'M5': M5.str(10), 'TR': [x.str(5) for x in TR], 'scale': scale.str(30)}, open(LOG + 'z5_arb_coeffs.json', 'w'))
fac = [arb(factorial(n)) for n in range(K+1)]
def coeffs(t):
    pw = [arb(1)]
    for j in range(K): pw.append(pw[-1]*(-t))
    return [sum((a[m+2*j]*pw[j]*fac[m+2*j]/(fac[j]*fac[m]) for j in range((K-1-m)//2+1)), arb(0)) for m in range(K)]
def ev(b, w):
    s = acb(0)
    for x in reversed(b): s = s*w + x
    return s
def der(b): return [b[m+1]*(m+1) for m in range(len(b)-1)]
# (3) circle certificate on a time grid, boxes covering the circle and the time step at once
NC = 720; nst = 21
tgrid = [arb(-k)/10000 for k in range(nst)]
dt = arb('0.00005')                                        # each grid time covers [t - dt, t + dt]
ok = True; worst = []; maxH1 = 0.0; maxH2 = 0.0
for t in tgrid:
    tb = arb(t.mid(), dt.mid())                            # time ball
    b = coeffs(tb); bd = der(b); bdd = der(bd)
    mH = mHp = mR = None
    for j in range(NC):
        th = arb(2*j + 1)/NC; cen = r0*acb(th).exp_pi_i()
        half = float(0.6*3.1416/NC)*1.05
        box = acb(arb(cen.real.mid(), half), arb(cen.imag.mid(), half))
        hv = ev(b, box); hpv = ev(bd, box)
        hv = acb(hv.real + arb(0, TR[0].upper()), hv.imag + arb(0, TR[0].upper()))
        hpv = acb(hpv.real + arb(0, TR[1].upper()), hpv.imag + arb(0, TR[1].upper()))
        x1 = abs(hv).lower(); x2 = abs(hpv).lower()
        if t == tgrid[0]:
            maxH1 = max(maxH1, float(abs(hpv).upper())); maxH2 = max(maxH2, float(abs(ev(bdd, box)).upper()) + float(TR[2].upper()))
        # radial derivative: Re(conj(H) H' w/|w|) with w the box (w/|w| enclosed by box/r0)
        rad = (hv.conjugate()*hpv*box/r0).real.lower()
        mH = x1 if mH is None or x1 < mH else mH; mHp = x2 if mHp is None or x2 < mHp else mHp; mR = rad if mR is None or rad < mR else mR
    good = (mH > 0) and (mHp > 0) and (mR > 0)
    ok = ok and good; worst.append((float(t.mid()), float(mH), float(mHp), float(mR)))
print('CERT (Arb) on |w| = 0.6 for all t in [-2.05e-3, 5e-5]: no zero of H_t, no zero of H_t\', exit set empty:', ok,
      ' min |H| >=', min(x[1] for x in worst), ' min |H\'| >=', min(x[2] for x in worst), ' min Re conj(H)H_r >=', min(x[3] for x in worst), ' time', round(time.time()-T0, 1), flush=True)
# (4) collision window
# rigorous winding numbers at t = 0 on |w| = 0.6 (no zeros on the circle certified above): arg increments as balls
b0 = coeffs(arb(0)); bd0 = der(b0); NW = 720
def wind_rig(b):
    pts = [r0*acb(arb(2*j)/NW).exp_pi_i() for j in range(NW)]
    v = [ev(b, p) for p in pts]; tot = arb(0)
    for j in range(NW):
        q = v[(j+1) % NW]/v[j]
        ang = q.arg()
        if not (abs(ang) < arb(3)): return None
        tot += ang
    return tot/(2*arb.pi())
arcf = 2*3.14159266*0.6/NW
print('winding validity at t = 0 (arg change per arc <= arc*max|f\'|/min|f| < pi): H: %.4f, H\': %.4f rad' % (arcf*maxH1/worst[0][1], arcf*maxH2/worst[0][2]))
print('Arb winding numbers on |w| = 0.6 at t = 0: H ->', wind_rig(b0).str(5), "  H' ->", wind_rig(bd0).str(5), flush=True)
tc_num = arb(json.load(open(LOG + 'lehmer_tc.json'))['tc'])
eps = arb('1e-12')
t1 = tc_num + eps; t2 = tc_num - eps
b1 = coeffs(t1); b2 = coeffs(t2)
def Hreal(b, w):                                          # rigorous enclosure of H_t(w) for real w (ball), incl. truncation
    v = ev(b, acb(w)).real
    return v + arb(0, TR[0].upper())
wstar = arb('-1.4365816786e-6')
neg = Hreal(b1, wstar); pL = Hreal(b1, arb('-0.3')); pR = Hreal(b1, arb('0.3'))
print('t1 = t_c + 1e-12: H(w*) =', neg.str(5), ' H(-0.3) =', pL.str(5), ' H(0.3) =', pR.str(5), ' => two real zeros:', bool(neg < 0 and pL > 0 and pR > 0))
# positivity of H_{t2} on [-0.3, 0.3] by adaptive bisection
stack = [(arb('-0.6'), arb('0.6'))]; nint = 0; pos = True
while stack:
    lo, hi = stack.pop(); mid = (lo + hi)/2; radw = (hi - lo)/2
    v = Hreal(b2, arb(mid.mid(), radw.upper()))
    nint += 1
    if v > 0: continue
    if radw < arb('1e-14') or nint > 200000: pos = False; break
    stack.append((lo, mid)); stack.append((mid, hi))
print('t2 = t_c - 1e-12: H_t2 > 0 on [-0.6, 0.6]:', pos, ' (', nint, 'intervals)  => (with 2 zeros in |w|<0.6 by the circle certificate) the pair is non-real at t2')
print('CERT t_c in (', t2.str(15), ',', t1.str(15), ')  time', round(time.time()-T0, 1))
