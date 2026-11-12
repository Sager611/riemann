"""Joint optimisation over R with first-order P2: for each R the one-piece optimum (P1 = sinh, Q Ritz deg 59) and
the exact best P2 for that (P1, Q, R): kappa2(R) = 1 - log(c1 - L'M^-1 L)/R."""
import sys, time
import numpy as np
from tp import *
from onepiece import best
cases = {'uncond_KRZ': (5/54 - 1e-7, 25/298 - 1e-7), 'GRC_KRZ': (1/8 - 1e-7, 1/8 - 1e-7)}
Rs = [float(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14]
K = 8; B = GBasis(K)
for name, (xi1, xi2) in cases.items():
    print(f"== {name}: xi1={xi1:.9f} xi2={xi2:.9f}")
    for R in Rs:
        t0 = time.time()
        c, h, lam = best(R, xi1)
        P1 = P_sinh(lam); Q = QCheb.from_h(h)
        C1 = c1(P1, Q, R, xi1)
        L = c12_vec(P1, Q, R, xi1, xi2, B); M, _ = c2_mat(Q, R, xi2, B, 120, 24, 24)
        d = L @ np.linalg.solve(M, L)
        k1 = 1 - np.log(C1) / R; k2 = 1 - np.log(C1 - d) / R
        print(f"  R={R:6.3f} lam={lam:.5f} c1={C1:.6e} (ritz {c:.6e})  kappa1={k1:.10f}  delta={d:.4e}  delta/c1={d/C1:.3e}  gain={k2-k1:.3e}  kappa2={k2:.10f} [{time.time()-t0:.1f}s]", flush=True)
