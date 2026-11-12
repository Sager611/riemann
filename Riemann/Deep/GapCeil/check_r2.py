"""validate r-general code against GapLambda's r=2 float scan (logs/bm/bmscan_r2_n14.txt: kappa=3.25, d0..d8) and r=1 (bmscan_r1_n18_D15: kappa=3.2645 d0)."""
import numpy as np
import gcore as G
import gcore_r as GR
for r, kap, th, ref in [(2, 3.25, 0.2499, [0.839545, 0.884610, 0.917006, 0.940582, 0.958457, 0.972252, 0.983070, 0.991683, 0.998596]),
                        (1, 3.2645, 0.2499, [0.811054, 0.888552, 0.933327, 0.955965, 0.969430])]:
    Nb = 20
    M0, MS, MSS, B = GR.mats_r(kap * np.pi, th, r, Nb, 36)
    out = []
    for d in range(len(ref)):
        Ir = np.linalg.matrix_power(B.I, 2 * r)[:Nb, :d + 1]
        Q, _ = np.linalg.qr(Ir)
        m0, ms, mss = (Q.T @ A @ Q for A in (M0, MS, MSS))
        out.append(G.best_v(m0, ms, mss, lo=0.6, hi=2.2)[0] / kap ** 2)
    print('r=%d kappa=%.4f theta=%.4f:' % (r, kap, th), ' '.join('d%d:%.6f(%.6f)' % (d, o, rf) for d, (o, rf) in enumerate(zip(out, ref))))
