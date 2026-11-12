"""Perturbative two-piece gain at the one-piece optima of Degree2Mol (params.json, read-only):
for fixed (P1, Q, R) the best P2 gives c = c1 - L^T M^{-1} L (c12 = g.L, c2 = g^T M g, g = P2'')."""
import json, sys, time, os
from fractions import Fraction
import numpy as np
from tp import *
D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Degree2Mol', 'params.json')))
def load(name):
    p = D[name]
    return (float(Fraction(p['xi'])), float(Fraction(p['R'])), float(Fraction(p['lam'])),
            [float(x) for x in p['h']], float(p['kappa_float']))
cases = [('Bernard_uncond', [25 / 298, None]), ('RSC', [None]), ('AT21_uncond', [25 / 298, None])]
res = [(int(a), int(b), int(c)) for a, b, c in (x.split(',') for x in sys.argv[1:])] or [(120, 24, 24)]
for name, th2s in cases:
    xi, R, lam, h, kf = load(name)
    P1 = P_sinh(lam); Q = QCheb.from_h(h)
    C1 = c1(P1, Q, R, xi); k1 = 1 - np.log(C1) / R
    print(f"== {name}: xi1={xi:.9f} R={R:.6f} lam={lam:.6f} deg Q={2*len(h)-1}  c1={C1:.8f} kappa1={k1:.10f} (params {kf:.10f})  Q(1)={Q(1.0):.3e}")
    for th2 in th2s:
        th2 = xi if th2 is None else th2
        for rs in res:
            for K in (2, 4, 8, 12):
                t0 = time.time()
                B = GBasis(K)
                L = c12_vec(P1, Q, R, xi, th2, B)
                M, asym = c2_mat(Q, R, th2, B, *rs)
                ev = np.linalg.eigvalsh(M)
                g = -np.linalg.solve(M, L); dl = -(g @ L)
                k2 = 1 - np.log(C1 - dl) / R
                print(f"  th2={th2:.9f} res={rs} K={K:2d}: delta=L'M^-1 L={dl:.6e}  dkappa={k2-k1:.4e}  kappa2={k2:.10f}"
                      f"  |g|max={np.abs(g).max():.2e} eig(M) [{ev.min():.2e},{ev.max():.2e}] asym {asym:.1e} [{time.time()-t0:.1f}s]")
                sys.stdout.flush()
