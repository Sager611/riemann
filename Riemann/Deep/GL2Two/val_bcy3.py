"""Best P2 in BCY's own ansatz span{x^3, x^4, x^5} (= GBasis(3)) for both BCY parameter sets; compare with the printed P2."""
import numpy as np
from tp import *
for name, Qz, P1c, P2c, R in [("kappa  (l.273-284)", [0.492, 0.604, 0, -0.08, 0, -0.06, 0, 0.046], [0, .842706, .00845721, .093117, .118788, -.0630687], [0, 0, 0, .0245412, -.00635566, .00603128], 1.28),
                              ("kappa* (l.285-294)", [0.485, 0.515], [0, .829473, .0104358, .082009, .177482, -.0993997], [0, 0, 0, .0323061, -.00553783, .0769594], 1.12)]:
    Qz = np.asarray(Qz, float); Qz = Qz / Qz.sum()
    Q = QCheb.from_z_mono(Qz); P1 = PMono(P1c); B = GBasis(3)
    C1 = c1(P1, Q, R, 4 / 7); L = c12_vec(P1, Q, R, 4 / 7, 1 / 2, B); M, _ = c2_mat(Q, R, 1 / 2, B, 100, 32, 32)
    g = -np.linalg.solve(M, L); P2 = B.P2_coeffs(g)
    k = lambda gg: 1 - np.log(C1 + 2 * gg @ L + gg @ M @ gg) / R
    gp = gvec_from_P2(P2c, 3)
    alt = list(P2c); alt[5] = -alt[5] / 10   # hypothesis: x^5 coefficient misprinted
    print(f"{name}: optimal P2 = {P2[3]:.7f} x^3 + {P2[4]:.8f} x^4 + {P2[5]:.7f} x^5 -> kappa {k(g):.7f};  printed P2 -> {k(gp):.7f};"
          f"  printed with x^5 coeff -> {alt[5]:.8f}: {k(gvec_from_P2(alt, 3)):.7f};  one piece {1-np.log(C1)/R:.7f}")
