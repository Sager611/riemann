"""Driver for gcert_r (general r, weighted Jacobi basis): certify that BM's criterion with this r fails at kappa for every
amplifier, every v and every theta <= 1/4. usage: runcert_r.py kappa N r vlo vhi [cell_shrink=0.7] [D]"""
import sys, time, json, resource
import numpy as np
from scipy.linalg import eigh
from flint import arb, fmpq
import gcert as GC
import gcert_r as GR

kap = GC.dec(sys.argv[1]); N = int(sys.argv[2]); RR = int(sys.argv[3]); vlo = GC.dec(sys.argv[4]); vhi = GC.dec(sys.argv[5])
SHRINK = float(sys.argv[6]) if len(sys.argv) > 6 else 0.7
t0 = time.time()
gam = int(float(2 * GC.A_(kap) * arb.pi().mid() + 3)) + 1
D = int(sys.argv[7]) if len(sys.argv) > 7 else max(50, int(2.5 * gam / 4 + 34))   # r >= 2: the tail bound carries ((i+2r)(j+2r))^(2r), so large kappa needs larger D
import math
from flint import ctx
GC.set_prec(128 + int(D * math.log2(max(gam / 4, 2))) + 32)   # Horner radius growth |I|^i |Gam_i| needs ~D log2(gam th) extra bits
ops = GR.OpsR(kap, N, D, gam, RR)
k2 = GC.A_(kap) ** 2; k2f = float(kap) ** 2; A = GC.A_
print('kappa=%s N=%d D=%d gam=%d Nw=%d prec=%d (Gamma %.0fs, setup %.0fs)' % (kap, N, D, gam, ops.Nw, ctx.prec, ops.t_gam, time.time() - t0), flush=True)
def npm(h, p):
    return np.array([[float(x.mid()) for x in row[:N]] for row in ops.A(h, p).tolist()[:N]])
F = {k: npm(*k) for k in [('a', 0), ('a', 1), ('a', 2), ('b', 0), ('b', 1), ('c', 0)]}
print('  matrices assembled (%.0fs)' % (time.time() - t0), flush=True)
M0 = F[('a', 0)]; MS = F[('a', 1)] + F[('b', 0)]; MSS = F[('a', 2)] + 2 * F[('b', 1)] + F[('c', 0)]
Gd = np.eye(N)
def top(v):
    return eigh(M0, v * v * M0 - 2 * v * MS + MSS, eigvals_only=True)[-1]
from scipy.optimize import minimize_scalar
r = minimize_scalar(lambda v: -top(v), bounds=(0.8, 2.0), method='bounded', options={'xatol': 1e-9})
vstar = r.x; lam = -r.fun
fr = eigh(MS, M0, eigvals_only=True)
print('  float: max_v lam = %.7f at v* = %.6f ; lam/k^2 = %.7f ; spec(M0^-1 MS) in [%.4f, %.4f]' % (lam, vstar, lam / k2f, fr[0], fr[-1]), flush=True)
res = {'kappa': str(kap), 'r': RR, 'N': N, 'D': D, 'gam': gam, 'prec': ctx.prec, 'float_lam_over_k2': lam / k2f, 'vstar': vstar, 'vlo': str(vlo), 'vhi': str(vhi), 'certs': []}
allok = True
def run(bl, label):
    global allok
    ok, info = GR.certify_r(ops, bl[(0, 0)], label)
    allok &= ok
    res['certs'].append({'ok': ok, **info})
    print('  %-44s %s deltaD=%s C=%s deltaW=%s chol=%s (%.0fs)' % (label, ok, info['deltaD'], info['C'], info['deltaW'], info.get('cholesky'), time.time() - t0), flush=True)
run({(0, 0): [('a', 0, arb(1))]}, 'M0 >= 0')
run({(0, 0): [('a', 1, arb(1)), ('b', 0, arb(1)), ('a', 0, -A(vlo))]}, 'MS - vlo M0 >= 0')
run({(0, 0): [('a', 1, arb(-1)), ('b', 0, arb(-1)), ('a', 0, A(vhi))]}, 'vhi M0 - MS >= 0')
def cell_min(vm, h):
    M1 = vm * vm * M0 - 2 * vm * MS + MSS; L = k2f * M1 - M0; Sp = MS - vm * M0
    return min(eigh(L + 2 * k2f * h * Sp, Gd, eigvals_only=True)[0], eigh(L - 2 * k2f * h * Sp, Gd, eigvals_only=True)[0])
def maxh(vm):
    lo, hi = 0.0, 0.8
    if cell_min(vm, 1e-9) <= 0: return 0.0
    for _ in range(36):
        mid = (lo + hi) / 2
        if cell_min(vm, mid) > 0: lo = mid
        else: hi = mid
    return lo
left = vlo; cells = []
while left < vhi and len(cells) < 60:
    lo, hi = 0.0, 0.8
    lf = float(left)
    for _ in range(28):
        mid = (lo + hi) / 2
        if maxh(lf + mid) >= mid: lo = mid
        else: hi = mid
    h = GC.dec('%.5f' % (lo * SHRINK))
    if h <= 0:
        print('  cell search FAILED at v=%s' % left); allok = False; break
    vm = left + h; cells.append((vm, h)); left = vm + h
print('  %d cells for [%s, %s]' % (len(cells), float(vlo), float(vhi)), flush=True)
for (vm, h) in cells:
    for sg in (1, -1):
        bl = {(0, 0): [('a', 0, k2 * A(vm) ** 2 - 1 - 2 * sg * k2 * A(h) * A(vm)), ('a', 1, -2 * k2 * A(vm) + 2 * sg * k2 * A(h)),
                       ('b', 0, -2 * k2 * A(vm) + 2 * sg * k2 * A(h)), ('a', 2, k2), ('b', 1, 2 * k2), ('c', 0, k2)]}
        run(bl, 'cell [%.5f,%.5f] %+d' % (float(vm - h), float(vm + h), sg))
covered = bool(cells) and cells[0][0] - cells[0][1] <= vlo and cells[-1][0] + cells[-1][1] >= vhi
res['cells'] = [[str(a), str(b)] for a, b in cells]
res['covered'] = covered; res['ALL_OK'] = bool(allok and covered)
res['peak_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20; res['seconds'] = time.time() - t0
print('RESULT kappa=%s: ALL_OK=%s covered=%s ; peak %.0f MB, %.0fs' % (kap, res['ALL_OK'], covered, res['peak_MB'], res['seconds']))
json.dump(res, open('logs/certR_r%d_%s.json' % (RR, sys.argv[1]), 'w'), indent=1)
print('JSON logs/certR_r%d_%s.json' % (RR, sys.argv[1]))
