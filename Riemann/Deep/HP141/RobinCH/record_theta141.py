# record_theta141.py -- RobinCH record claims (check_all2.py part (3), as STCAudit/dep_robinch_record.py) with PR 141's exact theta141,
# in this folder's copies (certify.py C3 = 9.4926, STCAudit's correction). Claims: 3-digit (as printed) and 4-digit roundings DOWN of the
# endpoints of logs/endpoint_level4_theta141.log (10^48.252564, 10^47.156287) and of logs/endpoint_level4_old.log (10^48.252447, 10^47.156206).
import sys, os, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from flint import arb, ctx
ctx.prec = 192
import certify3 as Q
import certify as C
import theta141 as T
print('certify.C3 =', C.C3, flush=True)
t0 = time.time(); LN10 = Q.LN10
e1 = Q.eta1_ch()
TH = f'{T.TH_NUM}/{T.TH_DEN}'
claims = [('theta141, 4-digit', TH, '1.788e48', '1.433e47'), ('old 0.87338196661320752, 4-digit', '0.87338196661320752', '1.788e48', '1.432e47'),
          ('theta141, 3-digit (printed record)', TH, '1.78e48', '1.43e47')]
for label, b0s, XN, XR in claims:
    b0 = Q.parse_b0(b0s)
    okN, LN, nbN, mN, kN = Q.run3('N', b0, 7*LN10, arb(XN).log(), e1, 4)
    okR, LR, nbR, mR, kR = Q.run3('R', b0, 10*LN10, arb(XR).log(), e1, 4)
    LX = arb(XR).log()
    eta = Q.best_eta3(LX, e1 if LX + arb('1e-3') > Q.LY1 else arb(0))
    log10n = float((arb(XR)*(1 - eta.upper())/LN10).lower())
    print(f'Level 3 (code 4), b0 = {label}: Nicolas on [1e7, {XN}]: {"PASS" if okN else "FAIL"} ({nbN} boxes, min {mN:.3g}, k <= {kN}); '
          f'Robin on CA with largest prime in [1e10, {XR}]: {"PASS" if okR else "FAIL"} ({nbR} boxes, min {mR:.3g}, k <= {kR}) '
          f'=> Robin for 5041 <= n <= 10^(10^{math.floor(math.log10(log10n)*10000)/10000:.4f}); {time.time()-t0:.0f}s', flush=True)
