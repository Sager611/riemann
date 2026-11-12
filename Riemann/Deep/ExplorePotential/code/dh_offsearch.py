# dh_offsearch.py (explore35): locate DH off-line zeros (sigma > 1/2) below TMAX from the on-line count deficit, close Landau's sum rule.
# Seeds: every gap between consecutive on-line zeros whose deficit increment exceeds THR; vectorized Newton on f (dh_ft.f), mpmath polish.
import sys, json, math, time, numpy as np, mpmath as mp
sys.path.insert(0, '.')
from dh_ft import f, vartheta, coeffs
d = json.load(open(sys.argv[1])); THR = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6
on = np.array(d['online']); T = d['TMAX']; M = int(T/2) + 60; t0 = time.time()
vt = vartheta(on); n = np.arange(1, len(on) + 1); dev = vt/np.pi - (n - 0.5)
inc = np.diff(dev); gi = np.nonzero(inc > THR)[0]
print('gaps', len(on) - 1, 'seeded', len(gi), 'total deficit %.2f' % dev[-1], flush=True)
lo = on[gi]; hi = on[gi + 1]
seeds = []
for sg in (0.53, 0.58, 0.66, 0.78, 0.95, 1.15):
    for fr in (0.2, 0.5, 0.8):
        seeds.append(sg + 1j*(lo + fr*(hi - lo)))
s = np.concatenate(seeds); s0 = s.copy()
h = 1e-6
for it in range(50):
    fs = f(s, M); dd = (f(s + h, M) - f(s - h, M))/(2*h)
    step = fs/dd; step = np.where(np.abs(step) > 0.3, 0.3*step/np.abs(step), step)
    s = s - step
    s = np.where((s.real < 0.45) | (s.real > 2.5) | (s.imag < 1) | (s.imag > T + 2), np.nan + 0j, s)
fin = np.isfinite(s)
val = np.abs(f(np.nan_to_num(s), M))
ok = fin & (val < 1e-8) & (s.real > 0.5 + 1e-6) & (s.imag <= T)
cand = s[ok]; cand = cand[np.argsort(cand.imag)]
uniq = []
for z in cand:
    if not uniq or min(abs(z - u) for u in uniq[-30:]) > 1e-5: uniq.append(z)
cand = np.array(uniq)
print('newton done', round(time.time() - t0, 1), 'converged off-line candidates', ok.sum(), 'unique', len(cand), flush=True)
mp.mp.dps = 25
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
def fmp(z): return mp.power(5, -z)*(mp.zeta(z, mp.mpf(1)/5) + kap*mp.zeta(z, mp.mpf(2)/5) - kap*mp.zeta(z, mp.mpf(3)/5) - mp.zeta(z, mp.mpf(4)/5))
roots = []
for r0 in cand:
    try: r = mp.findroot(fmp, mp.mpc(r0.real, r0.imag), tol=1e-18)
    except Exception: continue
    if float(r.real) > 0.5 + 1e-6 and all(abs(r - q) > 1e-8 for q in roots): roots.append(r)
roots = sorted(roots, key=lambda r: float(r.imag))
off = [[float(r.real), float(r.imag)] for r in roots]
print('distinct off-line right members', len(off), '-> zeros off the line', 2*len(off), 'vs deficit %.1f' % dev[-1], round(time.time() - t0, 1), flush=True)
bet = np.array([b for b, g in off]); gam = np.array([g for b, g in off])
print('delta = beta - 1/2: quartiles', np.round(np.quantile(bet - 0.5, [0, 0.25, 0.5, 0.75, 1]), 4).tolist(), 'count beta > 1:', int((bet > 1).sum()))
a, c = coeffs(40)
print('| x | Landau P | Re S_on - P | Re (S_on + S_off) - P | on-line share Re S_on/P | mean cos(gamma log x) over off-line |')
res = []
for x in range(2, 41):
    L = math.log(x); S_on = np.exp(1j*on*L).sum(); S_off = (np.exp(1j*gam*L)*2*np.cosh((bet - 0.5)*L)).sum()
    P = -(T/(2*math.pi))*c[x]/math.sqrt(x); res.append(abs((S_on + S_off).real - P))
    print('| %d | %+.1f | %+.1f | %+.1f | %s | %+.3f |' % (x, P, S_on.real - P, (S_on + S_off).real - P, ('%.2f' % (S_on.real/P)) if abs(P) > 1e-9 else '-', np.cos(gam*L).mean()))
print('mean |sum-rule residual| with off-line zeros: %.2f' % np.mean(res))
# phase locking histogram for x = 2, 3
for x in (2, 3, 4):
    ph = np.mod(gam*math.log(x), 2*math.pi); hist = np.histogram(ph, bins=8, range=(0, 2*math.pi))[0]
    print('x = %d: off-line phase gamma log x mod 2pi, 8 bins from 0:' % x, hist.tolist(), ' on-line:', np.histogram(np.mod(on*math.log(x), 2*math.pi), bins=8, range=(0, 2*math.pi))[0].tolist())
json.dump({'TMAX': T, 'offline': off}, open(sys.argv[2], 'w'))
print('done', round(time.time() - t0, 1))
