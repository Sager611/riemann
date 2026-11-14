# RobinCH record claims only (from check_all2.py parts (2)/(3)), with the patched certify.py (9.4926) in this folder.
import sys, time, math
sys.path.insert(0, '/tmp/stc2/RobinCH')
from flint import arb, ctx
ctx.prec = 192
import certify3 as Q
import certify as C
print('certify.C3 =', C.C3, flush=True)
t0 = time.time(); LN10 = Q.LN10
e1 = Q.eta1_ch()
claims = [('(3) Level 3', 4, '0.87338196661320752', '1.78e48', '1.43e47'), ('(2) Level 2-CH', 2, '0.87338196661320752', '7.17e41', '8.47e40')]
for part, level, label, XN, XR in claims:
    b0 = Q.parse_b0(label)
    okN, LN, nbN, mN, kN = Q.run3('N', b0, 7*LN10, arb(XN).log(), e1, level)
    okR, LR, nbR, mR, kR = Q.run3('R', b0, 10*LN10, arb(XR).log(), e1, level)
    LX = arb(XR).log()
    eta = Q.best_eta3(LX, e1 if LX + arb('1e-3') > Q.LY1 else arb(0))
    log10n = float((arb(XR)*(1 - eta.upper())/LN10).lower())
    print(f'{part}, b0 = {label}: Nicolas on [1e7, {XN}]: {"PASS" if okN else "FAIL"} ({nbN} boxes, min {mN:.3g}, k <= {kN}); '
          f'Robin on CA with largest prime in [1e10, {XR}]: {"PASS" if okR else "FAIL"} ({nbR} boxes, min {mR:.3g}, k <= {kR}) '
          f'=> Robin for 5041 <= n <= 10^(10^{math.floor(math.log10(log10n)*10000)/10000:.4f}); {time.time()-t0:.0f}s', flush=True)
