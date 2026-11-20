# hdx_nested.py (explore64): smallest nested cases (float64, as hdx_garland.py): links almost the window.
# N13: window lambda'^2 = 40, 41.5 from links lambda^2 = 227/6; N1: window 49, 50 from links 47. Hard and cosine translation links.
import sys, math, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreFunctional/code'); sys.path.insert(0, 'Riemann/Deep/ExploreHDX')
from xform import Form
from hdx_garland import Rmat
n = int(sys.argv[1]) if len(sys.argv) > 1 else 24
print(f'# n={n}/sector: c, c_link, sector, lambda_min R (hard), lambda_min R (cosine), lambda_min R on modes >= 8 (hard)')
for c, cl in [(40.0, 227/6), (41.5, 227/6), (49.0, 47.0), (50.0, 47.0)]:
    F = Form(c); ap = math.log(cl)/2
    for sec in ['even', 'odd']:
        Rh = Rmat(F, ap, n, sec, 'hard'); Rs = Rmat(F, ap, n, sec, 'smooth')
        lh = np.linalg.eigvalsh(Rh)[0]; ls = np.linalg.eigvalsh(Rs)[0]; hf = np.linalg.eigvalsh(Rh[8:, 8:])[0]
        print(f'{c:5.1f} {cl:8.4f} {sec:4s}  Rh {lh: .4f}  Rs {ls: .4f}  hf8 {hf: .4f}', flush=True)
