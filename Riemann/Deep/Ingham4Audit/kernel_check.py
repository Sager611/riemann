# Task (1): independent 256-bit check of the trapezoid kernel (DERIVATION sec. 3') and of TrapKernel's stored moments.
# G_a written from the definition (S(y) = max(sin y/y, s1) for y <= pi, min(s1, 1/y) for y >= pi), piecewise mpmath quadrature
# between the exact breakpoints, exact tail C2 lam^{p-2} p/(2-p) beyond lam = 93 (there G_a = C2/lam^2 exactly). cert.py is used only
# to build its TrapKernel (the object under test).
import os, sys, random, math
from mpmath import mp, mpf, sin, pi, quad, findroot, exp, log, mpc, fabs
mp.prec = 256
s1 = mpf('0.2173'); a = mpf(1)/10; b = 1 - a
def S(y):
    if y == 0: return mpf(1)
    return max(sin(y)/y, s1) if y <= pi else min(s1, 1/y)
def G(l, dm):
    if l == 0: return (1 + 2*dm*a)*(1 + 2*dm*b)
    return (S(a*l/2) + 2*dm*min(a, 1/l))*(S(b*l/2) + 2*dm*min(b, 1/l))
ys = findroot(lambda y: sin(y)/y - s1, mpf('2.6'))        # sinc = s1 on (0, pi)
yt = 1/s1                                                  # 1/y = s1
BRK = sorted(set([2*ys/a, 2*ys/b, 2*pi/a, 2*pi/b, 2*yt/a, 2*yt/b, 1/a, 1/b]))
YEND = mpf(93)
assert BRK[-1] < YEND
def C2(dm): return (2/a + 2*dm)*(2/b + 2*dm)
def Ddef(p, l0, dm):      # p int_{l0}^oo G lam^{p-1} dlam
    pts = [l0] + [x for x in BRK if x > l0] + ([YEND] if l0 < YEND else [])
    tot = mpf(0)
    for u, v in zip(pts, pts[1:]):
        tot += quad(lambda l: G(l, dm)*l**(p - 1), [u, v])
    top = max(l0, YEND)
    tot = p*tot + p*C2(dm)*top**(p - 2)/(2 - p)
    return tot
def Mom(p, l0, dm):       # max(1,l0)^p G(l0) + p int_{max(1,l0)}^oo G lam^{p-1}
    m = max(mpf(1), l0)
    return m**p*G(l0, dm) + Ddef(p, m, dm)
# kernel envelope at random z (both estimates): |Phi~(z)|/(ab) <= G(|Im z|)
random.seed(20261120)
worst = mpf(0)
for dm in (mpf('1e-9'), mpf('0.025'), mpf('0.15')):
    for t in range(3000):
        s = mpf(random.random())*dm; lam = mpf(10)**(mpf(random.uniform(-3, 4)))
        z = mpc(s, lam)
        ph = ((exp(a*z) - 1)/z)*((exp(b*z) - 1)/z)
        worst = max(worst, fabs(ph)/(a*b)/G(lam, dm))
print('envelope: max |Phi~_a(z)|/(ab G_a(Im z)) over 9000 random z (3 dmax, lam in [1e-3, 1e4]) = %s (must be <= 1)' % mp.nstr(worst, 12))
# G_a non-increasing and <= C2/lam^2 beyond 2 pi/a on a fine grid
mono = True; tail = mpf(0); prev = None
for j in range(0, 20001):
    l = mpf(j)/100
    for dm in (mpf('1e-9'), mpf('0.025')):
        pass
    v = G(l, mpf('0.025'))
    if prev is not None and v > prev: mono = False
    prev = v
    if l >= 2*pi/a: tail = max(tail, v*l**2/C2(mpf('0.025')))
print('G_a (dmax 0.025) non-increasing on lam = 0, 0.01, ..., 200: %s; max G_a lam^2/C2 for lam >= 2pi/a: %s' % (mono, mp.nstr(tail, 15)))
# TrapKernel under test
sys.dont_write_bytecode = True
import importlib.util
spec = importlib.util.spec_from_file_location('c4', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Ingham4', 'cert.py'))   # read-only import, no bytecode written
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
PG = [C.arb(i)/200 for i in range(1, 301)]
out = []
for dms in ('1e-9', '0.025'):
    K = C.TrapKernel(C.arb('0.1'), C.arb(dms), PG)
    dm = mpf(dms)
    print('TrapKernel dmax=%s: N=%d, Lambda=%s, istore=%d, lam_istore=%s, C2=%s, G0=%s' % (dms, K.N, K.lam[K.N].str(8, radius=False), K.istore, K.lam[K.istore].str(8, radius=False), K.C2.str(10, radius=False), K.G0.str(10, radius=False)))
    # 5 stored entries per kernel: (j, i)
    for j, i in [(1, 0), (80, 555), (120, 300), (208, 1200), (300, K.istore)] if dms == '0.025' else [(201, 0), (72, 555), (100, 2000), (250, 40), (60, K.istore)]:
        p = mpf(j)/200; li = mpf(1.002)**0  # placeholder
        li = exp(i*log(1 + mpf(1)/500))
        d_ = Ddef(p, li, dm); st = mpf(K.T[str(PG[j - 1])][i])
        gstored = mpf(C.up(K.Gv[i]).mid().str(60, radius=False)); gdef = G(li, dm)
        out.append((dms, j, i, st/d_, gstored >= gdef))
        print('  T[p=%d/200][i=%d] (lam_i=%s): stored %s, 256-bit definition %s, stored/def = %s; G stored >= G def: %s' % (j, i, mp.nstr(li, 8), mp.nstr(st, 12), mp.nstr(d_, 12), mp.nstr(st/d_, 8), gstored >= gdef))
    # the three mom() branches against the definition of Mom
    for j, l0s in [(80, '0.5'), (80, '3.0294'), (120, '64.5'), (300, '100'), (60, '1000')]:
        p = mpf(j)/200; l0 = mpf(l0s)
        mc = mpf(K.mom(PG[j - 1], C.arb(l0s)).mid().str(60, radius=False)); md = Mom(p, l0, dm)
        print('  mom(p=%d/200, l0=%s): cert %s, definition %s, ratio %s' % (j, l0s, mp.nstr(mc, 12), mp.nstr(md, 12), mp.nstr(mc/md, 8)))
        out.append((dms, j, l0s, mc/md, mc >= md))
print('all stored/def >= 1 and all mom ratios >= 1:', all(o[3] >= 1 for o in out) and all(o[4] for o in out))
