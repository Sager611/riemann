# Arb (python-flint) checks: N(T0); Lemma 3.3; Prop 5.1 extremal cell; Pintz Lemma 4 sanity with Gallego's zeros.
import sys, time, math, os
from flint import arb, acb, ctx
DATA = sys.argv[1]
t0 = time.time()
# 1. N(T0)
ctx.prec = 64
N = arb(1300000).zeta_nzeros(); print('1. Arb N(1.3e6) =', N, ' (%.1fs)' % (time.time() - t0), flush=True)
# 2. Lemma 3.3: min |zeta(1/4+it)| on [0,4], 2^13 balls of radius 2^-12 (Arb lower bounds)
t1 = time.time(); m = None; h = arb(4) / 8192
for j in range(8192):
    tb = arb(h * j + h / 2, h / 2)          # ball covering [jh, (j+1)h]
    z = acb(arb(1) / 4, tb).zeta()
    lo = abs(z).lower()
    if m is None or lo < m: m = lo; jm = j
print('2. min_j lower|zeta(1/4+it)| over 8192 balls =', m.str(8), 'at t ~', float(h.mid()) * jm, ' (%.1fs)' % (time.time() - t1), flush=True)
# 3. Prop 5.1 extremal cell nu = 5042995 (Bi-Hou Remark 5.2): zeros nu-6..nu+6 and N_nu for nu-2..nu+2
t1 = time.time(); ctx.prec = 128
nu = 5042995; G = {}
for n in range(nu - 6, nu + 7):
    G[n] = acb.zeta_zero(n).imag
print('3. zeros', nu - 6, '..', nu + 6, 'computed (%.1fs)' % (time.time() - t1), flush=True)
for n in (nu - 1, nu, nu + 1):
    print('   gamma_%d = %s' % (n, G[n].str(15)))
def cellN(v):
    lo = (G[v - 1] + G[v]) / 2; hi = (G[v] + G[v + 1]) / 2; s = arb(0)
    for mu in range(nu - 6, nu + 7):
        if mu == v: continue
        d = lo - G[mu] if mu < v else G[mu] - hi
        if d < 1: s += d.log()          # dist(gamma_mu, I_v) < 1
    return lo, hi, s
for v in range(nu - 2, nu + 3):
    lo, hi, s = cellN(v)
    print('   cell %d: [%s, %s]  N_nu = %s' % (v, lo.str(12), hi.str(12), s.str(10)))
print('   (Bi-Hou: N_nu >= -8.76707 overall, extremal at nu=5042995; theta_{nu-1}~2650984.27884, theta_nu~2650984.61215)', flush=True)
# 4. Pintz Lemma 4 sanity: |log|zeta(s)| - sum_{|gamma-t|<c0} log|s-rho|| vs 11 log t + 14, zeros from Gallego's data
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crho_format import CrhoFile
c0 = math.tan(1.0); ctx.prec = 64
cf = CrhoFile(os.path.join(DATA, 'cc_00742987.crho'))   # n = 742987..1114479, gamma ~ 4.6e5..6.6e5
gam = [float(cf.txt('gamma', i)) for i in range(0, cf.count)]
del cf
import bisect
worst = 0.0
for t in (500000.123, 550000.0, 600000.777, 650000.31):
    for sig in (0.5 + 1e-4, 0.5 + 1e-2, 0.75, 1.0, 1.5, 2.0):
        z = acb(sig, t).zeta(); lz = float(abs(z).log().mid())
        i0 = bisect.bisect_left(gam, t - c0); i1 = bisect.bisect_right(gam, t + c0)
        sl = sum(0.5 * math.log((sig - 0.5) ** 2 + (t - g) ** 2) for g in gam[i0:i1])
        d = abs(lz - sl); worst = max(worst, d / (11 * math.log(t) + 14))
        print('   t=%.3f sigma=%.4f zeros=%d  |log|z|-sum| = %.3f  vs 11log t+14 = %.1f' % (t, sig, i1 - i0, d, 11 * math.log(t) + 14))
print('4. worst ratio to Pintz bound:', round(worst, 4), ' total %.1fs' % (time.time() - t0))
