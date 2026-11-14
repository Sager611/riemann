# RobinHeight record claims only (from check_all.py parts (2) and (4)), with the patched certify.py (9.4926) in this folder.
import sys, time, math
sys.path.insert(0, '/tmp/stc2/RobinHeight')
from flint import arb, ctx
ctx.prec = 192
import certify as C
import certify2 as K
print('certify.C3 =', C.C3, flush=True)
t0 = time.time(); LN10 = arb(10).log()
for label, b0, XN, XR in [('7/8', arb(7)/8, '1.04e40', '9.37e38'), ('0.87338196661320752', arb('0.87338196661320752'), '1.31e40', '1.19e39')]:
    e1 = K.eta1_of(b0)
    okN, LN, nbN, mN = K.run2('N', b0, 7*LN10, arb(XN).log(), e1)
    okR, LR, nbR, mR = K.run2('R', b0, 10*LN10, arb(XR).log(), e1)
    LX = arb(XR).log()
    eta = K.best_eta(LX, b0, e1 if LX > K.LY1 else arb(0))
    log10n = float((arb(XR)*(1 - eta.upper())/LN10).lower())
    print(f'(4) Level 2, b0 = {label}: Nicolas on [1e7, {XN}]: {"PASS" if okN else "FAIL"} ({nbN} boxes, min {mN:.3g}); '
          f'Robin on CA with largest prime in [1e10, {XR}]: {"PASS" if okR else "FAIL"} ({nbR} boxes, min {mR:.3g}); '
          f'=> Robin for 5041 <= n <= 10^(10^{math.log10(log10n):.4f}); {time.time()-t0:.0f}s', flush=True)
for label, b0, dlt, XN, XR in [('7/8', arb(7)/8, arb(10)**-9, '1.93e32', '1.56e31'), ('0.87338196661320752', arb('0.87338196661320752'), arb(10)**-9, '2.67e32', '2.13e31')]:
    C.setdelta(dlt)
    okN, LN, nbN, mN = C.run('N', b0, 7*LN10, arb(XN).log())
    okR, LR, nbR, mR = C.run('R', b0, 10*LN10, arb(XR).log())
    print(f'(2) Level 1, b0 = {label}: Nicolas on [1e7, {XN}]: {"PASS" if okN else "FAIL"} ({nbN} boxes, min {mN:.3g}); Robin [1e10, {XR}]: {"PASS" if okR else "FAIL"} ({nbR} boxes, min {mR:.3g}); {time.time()-t0:.0f}s', flush=True)
