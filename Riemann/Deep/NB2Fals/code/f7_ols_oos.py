# nb2fals f7: identifiable free-amplitude OLS fits (K = 0, 5; C2 + C3/L; weights 1/N) on N <= 5500 and their out-of-sample tests.
import os, sys, math, numpy as np
os.environ['NB2F_NMAX'] = '5500'
src = open('/Users/adriansager/playground/riemann/Riemann/Deep/NB2Fals/code/f1_fits.py').read().split("b20 = lambda m")[0].replace("NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4000", "NMAX = 5500")
exec(src)
for (lo, hi) in ((100, 640), (100, 5500)):
    m = sel(lo, hi); uu, NN, y = uall[m], Nall[m], Yall['C2'][m]
    for K in (0, 5):
        beta, a, res, cond, X = ols(uu, NN, y, K, True)
        se = [boot(uu, res, a, b) for b in (20 / np.median(NN), 0.6)]
        print(f'OLS C2 + C3/L K={K} [{lo},{hi}]: c = {beta[0]:+.5f}, C3 = {beta[1]:+.4f}; block-bootstrap s.e. {se[0]:.4f} (20-N blocks), {se[1]:.4f} (b = 0.6); cond {cond:.1e}', flush=True)
for (lo, mid, hi) in ((100, 400, 640), (100, 1280, 5500), (100, 4000, 5500)):
    m1 = sel(lo, mid); m2 = (Nall > mid) & (Nall <= hi)
    for K in (0, 5):
        beta, a, res, cond, X = ols(uall[m1], Nall[m1], Yall['C2'][m1], K, True)
        pred = design(uall[m2], K, True) @ beta; y2 = Yall['C2'][m2]; w2 = (1 / Nall[m2]) / np.sum(1 / Nall[m2])
        print(f'OOS OLS K={K}: fit [{lo},{mid}] -> ({mid},{hi}]: weighted mean error {w2 @ (y2 - pred):+.4f}, rms error {math.sqrt(w2 @ (y2 - pred) ** 2):.4f} vs sd(y) {math.sqrt(w2 @ (y2 - w2 @ y2) ** 2):.4f}; c(fit) = {beta[0]:+.4f}', flush=True)
