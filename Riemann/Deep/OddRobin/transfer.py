#!/usr/bin/env python3
# transfer.py -- odd Robin at odd CA numbers N_eps with x_1(eps) in [XD, XR] from the RobinCH Level-3 (level 4 code) Robin certificate.
#   uv run --no-project --with python-flint --with mpmath --with numpy python3 transfer.py [XD] [XR] [b0]
# For CA M = M_eps and odd CA N = N_eps = M / 2^a, a = alpha_2(eps):
#   odd Robin at N  <=>  R(M) := gamma + log log log M - log(sigma(M)/M)  >  C_odd := [log log log M - log log log N] - log(1 - 2^(-a-1)).
# RobinCH's margin m on a box of x = x_1 satisfies (sqrt x log x) R(M_eps) >= m for every eps with x_1(eps) in the box.
# Bounds (PROOF.md, Lemma O2): 2^(-a) < 2/x;  a log 2 <= log(((x+1) log x/log 2 + 2)/2);  log N >= theta(x) - log 2 >= x(1 - 1/log x) - log 2;
#   C_odd <= 2^(-a) + a log 2/(log N log log N).   So it suffices that m > (sqrt x log x) * [2/x + A/(LN log LN)] on the box.
import sys, os, time, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'RobinCH'))
from flint import arb, ctx
import certify3 as Q
C = Q.C
ctx.prec = 192
t0 = time.time()
LOG2 = arb(2).log()
def corr(La, Lb):
    L = C.box(La, Lb); x = L.exp()
    A = (((x + 1)*L/LOG2 + 2)/2).log()
    LN = x*(1 - 1/L) - LOG2
    return ((L/2).exp()*L*(2/x + A/(LN*LN.log()))).upper()
def run(b0, L0, L1, e1, w0=arb('0.25'), wmin=arb('1e-6'), verbose=False):
    La = L0; minm = None; nbox = 0; w = w0; kmax = 0
    while La < L1:
        Lb = La + w
        if Lb > L1: Lb = L1
        m, arg = Q.margins4('R', La, Lb, b0, e1); nbox += 1
        mo = m - corr(La, Lb)
        if verbose: print(f'  box [{float(La.mid()):.4f}, {float(Lb.mid()):.4f}] m = {float(m.lower()):.5g}, corr = {float(corr(La, Lb).mid()):.3g}', flush=True)
        if mo > 0:
            if minm is None or mo.lower() < minm: minm = float(mo.lower())
            if arg: kmax = max(kmax, arg[0])
            La = Lb
            if w < w0: w = w*2
        else:
            w = w/2
            if w < wmin: return False, La, nbox, minm, kmax
    return True, L1, nbox, minm, kmax
if __name__ == '__main__':
    XD = sys.argv[1] if len(sys.argv) > 1 else '5e8'
    XR = sys.argv[2] if len(sys.argv) > 2 else '1.43e47'
    b0 = Q.parse_b0(sys.argv[3] if len(sys.argv) > 3 else '0.87338196661320752')
    e1 = Q.eta1_ch()
    if len(sys.argv) > 4:   # probe mode: a few boxes
        for Ls in sys.argv[4].split(','):
            La = arb(Ls); Lb = La + arb('0.002')
            m, arg = Q.margins4('R', La, Lb, b0, e1)
            print(f'probe L = [{Ls}, +0.002]: m in [{float(m.lower()):.5g}, {float(m.upper()):.5g}], corr = {arb(corr(La, Lb)).str(4)}, k,lam = {arg}')
        sys.exit(0)
    L0 = arb(XD).log(); L1 = arb(XR).log()
    ok, Lend, nb, mn, km = run(b0, L0, L1, e1)
    LX = arb(XR).log()
    eta = Q.best_eta3(LX, e1 if LX + arb('1e-3') > Q.LY1 else arb(0))          # |theta(X_R)/X_R - 1| <= eta (RobinCH Lemma 12)
    lg10 = float(((arb(XR)*(1 - eta.upper()) - LOG2)/Q.LN10).lower())             # log N* >= theta(X_R) - log 2, N* = last odd CA checked
    Y = math.floor(math.log10(lg10)*10000)/10000
    print(f'{"PASS" if ok else "FAIL"}: odd Robin at every odd CA number with x_1 in [{XD}, {XR}] (b0 = {b0.str(10)}) => all odd N < n <= 10^(10^{Y:.4f}): '
          f'{nb} boxes, min margin after odd correction {mn:.3g}, k <= {km}' + ('' if ok else f'; stopped at x = e^{float(Lend.mid()):.4f}') + f'; {time.time()-t0:.0f}s')
    sys.exit(0 if ok else 1)

