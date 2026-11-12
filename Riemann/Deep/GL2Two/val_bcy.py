"""Validate tp.py against BCY 1002.4127 l.271-284 (kappa >= .4105 at th1=4/7, th2=1/2, R=1.28)
and KRZ 1609.03738 l.1750-1772 (kappa >= .0693872 at R_KRZ=2.82505, nu1=nu2=1/4; .0297607 at R_KRZ=3.21, nu=5/27, 25/149)."""
import numpy as np, sys, time
from tp import *

def run(name, Qz, P1c, P2c, R, th1, th2, K=6, res=(80, 24, 24), norm=None):
    t0 = time.time()
    Qz = np.asarray(Qz, float)
    if norm == 'scale': Qz = Qz / Qz.sum()                 # Q(0) = sum of (1-2x)^k coefficients
    if norm == 'shift': Qz = Qz.copy(); Qz[0] += 1 - Qz.sum()
    Q = QCheb.from_z_mono(Qz); P1 = PMono(P1c)
    name = name + (f' [Q(0)->1 by {norm}]' if norm else ' [as printed]')
    B = GBasis(K)
    g = gvec_from_P2(P2c, K)
    # check representation
    s = np.linspace(0, 1, 11)
    gs = np.polynomial.polynomial.polyval(s, np.polynomial.polynomial.polyder(np.asarray(P2c, float), 2))
    err = np.abs(g @ B.ev(s, None) - gs).max()
    C1 = c1(P1, Q, R, th1)
    L = c12_vec(P1, Q, R, th1, th2, B)
    M, asym = c2_mat(Q, R, th2, B, *res)
    c12 = g @ L; c2 = g @ M @ g
    c = C1 + 2 * c12 + c2
    k1 = 1 - np.log(C1) / R; k = 1 - np.log(c) / R
    gopt = -np.linalg.solve(M, L); copt = C1 + 2 * gopt @ L + gopt @ M @ gopt
    print(f"{name}: Q(0)={Q(0.0):.6f} P1(1)={P1(1.0):.7f} g-rep err {err:.1e} asym {asym:.1e}")
    print(f"  c1={C1:.10f} c12={c12:.6e} c2={c2:.6e} c={c:.10f}")
    print(f"  kappa one-piece {k1:.7f}  two-piece (their P2) {k:.7f}  (gain {k-k1:.3e});"
          f"  best P2 (K={K}) {1-np.log(copt)/R:.7f} (gain {1-np.log(copt)/R-k1:.3e})  [{time.time()-t0:.1f}s]")
    sys.stdout.flush()
    return C1, L, M

if __name__ == "__main__":
  for norm in (None, 'scale', 'shift'):
      # BCY
      run("BCY zeta 4/7,1/2", [0.492, 0.604, 0, -0.08, 0, -0.06, 0, 0.046],
          [0, .842706, .00845721, .093117, .118788, -.0630687],
          [0, 0, 0, .0245412, -.00635566, .00603128], 1.28, 4 / 7, 1 / 2, norm=norm)
      run("BCY zeta (res up)", [0.492, 0.604, 0, -0.08, 0, -0.06, 0, 0.046],
          [0, .842706, .00845721, .093117, .118788, -.0630687],
          [0, 0, 0, .0245412, -.00635566, .00603128], 1.28, 4 / 7, 1 / 2, res=(100, 32, 32), norm=norm)
      # KRZ, Ramanujan
      run("KRZ GRC nu=1/4,1/4", [.498939, 1.53685, 0, -2.7925, 0, 2.77524, 0, -1.01853],
          [0, .921756, .150879, -.371912, .488862, -.189585],
          [0, 0, 0, -.0000537029, .0000752763, -.000142568], 2 * 2.82505, 1 / 8, 1 / 8, norm=norm)
      # KRZ, unconditional
      run("KRZ uncond nu=5/27,25/149", [.499386, 1.58992, 0, -2.99061, 0, 3.01825, 0, -1.11694],
          [0, .93271, .147723, -.35572, .444208, -.168921],
          [0, 0, 0, -.0000665503, -.00016405, .0000736009], 2 * 3.21, 5 / 54, 25 / 298, norm=norm)
