"""Resolution check of the first-order gain at the one-piece optima (Degree2Mol params.json)."""
import json, sys, time, os
from fractions import Fraction
import numpy as np
from tp import *
D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Degree2Mol', 'params.json')))
for name, th2x in [('Bernard_uncond', 25/298 - 1e-7), ('RSC', None), ('AT21_uncond', 25/298 - 1e-7)]:
    p = D[name]; xi = float(Fraction(p['xi'])); R = float(Fraction(p['R'])); lam = float(Fraction(p['lam']))
    h = [float(s) for s in p['h']]; th2 = xi if th2x is None else th2x
    P1 = P_sinh(lam); Q = QCheb.from_h(h); C1 = c1(P1, Q, R, xi)
    for n12, res in [(40, (120, 24, 24)), (60, (160, 32, 32)), (80, (200, 40, 40))]:
        t0 = time.time(); B = GBasis(6)
        L = c12_vec(P1, Q, R, xi, th2, B, n=n12); M, _ = c2_mat(Q, R, th2, B, *res)
        g = -np.linalg.solve(M, L); d = -(g @ L)
        print(f"{name} th2={th2:.9f} n12={n12} res={res}: L[:3]={L[:3]}  M00={M[0,0]:.10e}  delta={d:.8e}  dkappa={np.log(C1/(C1-d))/R:.6e}  g[:3]={g[:3]} [{time.time()-t0:.0f}s]", flush=True)
