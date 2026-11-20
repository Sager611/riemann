# z7_dh_arb.py — explore66 adversarial test, rigorous (Arb): Davenport–Heilbronn's off-line pair under DH's heat flow.
# Same certificate as z5_arb.py (zeta): rigorous Taylor coefficients at c (fixed decimal centre on R), aliasing bound from
# max_{|w|=6}|H|, truncation majorant for |t| <= 0.2, |w| <= 1.2; circle |w| = 1.2: no zero of H_t, H_t', empty exit set for all
# t in [0, 0.2025] (time balls); winding 2 and 1 at t = 0; landing window t* +- 1e-12 (real pair after, H_t > 0 on R-diameter before).
import json, time, sys
from math import factorial
from flint import arb, acb, ctx, dirichlet_char
ctx.prec = 200
T0 = time.time()
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
chi = dirichlet_char(5, 2); chib = dirichlet_char(5, 3)
kap = ((arb(10) - 2*arb(5).sqrt()).sqrt() - 2)/(arb(5).sqrt() - 1)
c1 = (acb(1) - acb(0, 1)*kap)/2; c2 = c1.conjugate()
PI = acb(arb.pi()); I = acb(0, 1)
def H0(z):
    s = acb(0.5) + I*z/2
    return ((acb(5)/PI)**((s+1)/2))*((s+1)/2).gamma()*(c1*acb.dirichlet_l(s, chi) + c2*acb.dirichlet_l(s, chib))/8
c = arb('171.3986969707552')
R = arb(3); M = 128; K = M//2; rho = arb(6)
ws = [R*acb(arb(2*j)/M).exp_pi_i() for j in range(M)]
vals = [H0(c + w) for w in ws]
A = []
for k in range(K):
    s = acb(0)
    for j in range(M): s += vals[j]*acb(arb(-2*j*k)/M).exp_pi_i()
    A.append(s/(M*R**k))
NB = 3000; Mr = arb(0)
for j in range(NB):
    cen = rho*acb(arb(2*j + 1)/NB).exp_pi_i(); half = float(6*3.1416/NB)*1.05
    v = abs(H0(c + acb(arb(cen.real.mid(), half), arb(cen.imag.mid(), half)))).upper()
    if v > Mr: Mr = arb(v)
scale = A[0]
alias = Mr*(R/rho)**M/(1 - (R/rho)**M)
A = [(A[k] + acb(arb(0, (alias*rho**(-k)).upper()), arb(0, (alias*rho**(-k)).upper())))/scale for k in range(K)]
a = [x.real for x in A]
print('Mrho =', Mr.str(5), ' max |Im A_k| <=', max(float(abs(x.imag).upper()) for x in A), ' max rad a_k:', max(float(x.rad()) for x in a), ' time', round(time.time()-T0, 1), flush=True)
smax = arb('0.2025'); r0 = arb('1.2'); sig = (2*smax).sqrt(); sabs = abs(scale).lower()
# growth bound for ALL radii: for Re s >= 1/2, |L(s,chi)| <= sqrt2 |s|/Re s (partial summation, |sum_{n<=x} chi(n)| <= sqrt2), |Gamma(x+iy)| <= Gamma(x),
# and Lambda_f(s) = Lambda_f(1-s) (checked below). On |w| = rho around c: sigma' = max(Re s, 1 - Re s) <= 1/2 + rho/2, |s|, |1-s| <= 1/2 + (c + rho)/2.
s_t = acb('0.3', '50.0'); fe = abs(H0(-2*acb(0,1)*(s_t - acb(0.5))) ) ; fe2 = abs(H0(-2*acb(0,1)*((1 - s_t) - acb(0.5))))
print('functional equation check |Lambda(s)| vs |Lambda(1-s)| at s = 0.3+50i:', fe.str(10), fe2.str(10), flush=True)
cab = abs(c1) + abs(c2)
def G(rho):          # rigorous upper bound for max_{|w|=rho} |H0(c+w)|
    sp = arb(1)/2 + rho/2
    gam = ((sp + 1)/2).gamma(); gm = gam if gam > arb('1.2255') else arb('1.2255')
    return ((arb(5)/arb.pi())**((sp + 1)/2))*gm*cab*arb(2).sqrt()*2*(arb(1)/2 + (c + rho)/2)/8
import math
def best_rho(k):     # float search for a good radius (any radius gives a valid bound)
    best = None
    for i in range(0, 400):
        rho = 6*1.04**i
        sp = 0.5 + rho/2
        lg = ((sp + 1)/2)*math.log(5/math.pi) + max(math.lgamma((sp + 1)/2), math.log(1.2255)) + math.log(0.5 + (171.4 + rho)/2) - k*math.log(rho)
        if best is None or lg < best[0]: best = (lg, rho)
    return best[1]
def logterm(k, m):
    rho = arb(best_rho(k))
    return (G(rho)/sabs).log() - k*rho.log() + (arb(k+1).lgamma() - arb(k-m+1).lgamma()) + (k-m)*(r0 + sig*arb(k-m+1).sqrt()).log()
# also the rigorous radius-6 Cauchy bound (box cover Mr) for small k
def logterm6(k, m):
    return (Mr/sabs).log() - k*rho.log() + (arb(k+1).lgamma() - arb(k-m+1).lgamma()) + (k-m)*(r0 + sig*arb(k-m+1).sqrt()).log()
TR = []
for m in range(4):
    tot = arb(0); k = K; lastl = None
    while True:
        l1 = logterm(k, m); l2 = logterm6(k, m); l = l1 if l1 < l2 else l2
        tot += l.exp()
        if k > K + 50 and l < arb(-300):
            lnext = logterm(k+1, m); q = (lnext - l).exp()
            if q < arb('0.5'): break
        k += 1
        if k > 20000: break
    rats = [(logterm(kk+1, m) - logterm(kk, m)).exp().str(3) for kk in (k, 2*k, 4*k, 8*k)]
    far = (lnext.exp())/(1 - q)
    TR.append(tot + far)
    print('m=%d: tail <= %s (explicit to k=%d, remainder <= %s; successive ratios at k, 2k, 4k, 8k: %s)' % (m, (tot + far).str(4), k, far.str(3), rats), flush=True)
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
NC = int(sys.argv[1]) if len(sys.argv) > 1 else 1440; nst = int(sys.argv[2]) if len(sys.argv) > 2 else 41
dt = arb('0.2')/(nst - 1)/2
ok = True; mins = [9e9, 9e9, 9e9]
for i in range(nst):
    tb = arb((2*dt*i).mid(), dt.mid()*1.0000001)
    b = coeffs(tb); bd = der(b)
    for j in range(NC):
        cen = r0*acb(arb(2*j + 1)/NC).exp_pi_i(); half = float(1.2*3.1416/NC)*1.05
        box = acb(arb(cen.real.mid(), half), arb(cen.imag.mid(), half))
        hv = ev(b, box); hpv = ev(bd, box)
        hv = acb(hv.real + arb(0, TR[0].upper()), hv.imag + arb(0, TR[0].upper()))
        hpv = acb(hpv.real + arb(0, TR[1].upper()), hpv.imag + arb(0, TR[1].upper()))
        x1 = float(abs(hv).lower()); x2 = float(abs(hpv).lower()); x3 = float((hv.conjugate()*hpv*box/r0).real.lower())
        mins = [min(mins[0], x1), min(mins[1], x2), min(mins[2], x3)]
    ok = ok and min(mins) > 0
print('CERT (Arb) on |w| = 1.2 for all t in [-%s, 0.2+%s]: no zero of H_t, H_t\', exit set empty:' % (dt.str(3), dt.str(3)), ok, ' mins', mins, ' time', round(time.time()-T0, 1), flush=True)
b0 = coeffs(arb(0)); bd0 = der(b0); NW = 1440
def wind_rig(b):
    pts = [r0*acb(arb(2*j)/NW).exp_pi_i() for j in range(NW)]
    v = [ev(b, p) for p in pts]; tot = arb(0)
    for j in range(NW):
        ang = (v[(j+1) % NW]/v[j]).arg()
        if not (abs(ang) < arb(3)): return None
        tot += ang
    return tot/(2*arb.pi())
bdd0 = der(bd0)
B1u = sum((abs(bd0[mm])*r0**mm for mm in range(len(bd0))), arb(0)) + TR[1]; B2u = sum((abs(bdd0[mm])*r0**mm for mm in range(len(bdd0))), arb(0)) + TR[2]
arcw = 2*arb.pi()*r0/NW
print('winding validity (arg change per arc <= arc*max|f\'|/min|f|): H: %s, H\': %s rad' % ((arcw*B1u/arb(mins[0])).str(4), (arcw*B2u/arb(mins[1])).str(4)))
print('Arb winding on |w| = 1.2 at t = 0: H ->', wind_rig(b0).str(5), "  H' ->", wind_rig(bd0).str(5), flush=True)
tsn = arb(json.load(open(LOG + 'dh_tstar.json'))['tstar']); eps = arb('1e-12')
t1 = tsn + eps; t2 = tsn - eps; b1 = coeffs(t1); b2 = coeffs(t2)
def Hreal(b, w): return ev(b, acb(w)).real + arb(0, TR[0].upper())
wst = arb('-0.163668182664')
v1, vL, vR = Hreal(b1, wst), Hreal(b1, wst - arb('0.05')), Hreal(b1, wst + arb('0.05'))
print('t* + 1e-12: H(w*) =', v1.str(5), ' H(w*-0.05) =', vL.str(5), ' H(w*+0.05) =', vR.str(5), ' => real pair:', bool(v1 < 0 and vL > 0 and vR > 0))
# Taylor-shift the polynomial to the saddle w* so that interval Horner evaluation is second order near the minimum
from math import comb
bs = [sum((b2[mm]*arb(comb(mm, jj))*wst**(mm - jj) for mm in range(jj, len(b2))), arb(0)) for jj in range(len(b2))]
def Hsh(h): return ev(bs, acb(h)).real + arb(0, TR[0].upper())
stack = [(arb('-1.2') - wst, arb('1.2') - wst)]; nint = 0; pos = True
while stack:
    lo, hi = stack.pop(); mid = (lo + hi)/2; rw = (hi - lo)/2
    v = Hsh(arb(mid.mid(), rw.upper())); nint += 1
    if v > 0: continue
    if rw < arb('1e-14') or nint > 300000: pos = False; break
    stack.append((lo, mid)); stack.append((mid, hi))
print('t* - 1e-12: H > 0 on [-1.2, 1.2]:', pos, '(', nint, 'intervals) => conjugate pair; CERT t* in (', t2.str(15), ',', t1.str(15), ')  => Lambda_DH >', t2.str(12), ' time', round(time.time()-T0, 1))
