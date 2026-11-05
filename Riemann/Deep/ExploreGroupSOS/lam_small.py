import sys, math
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ChekAudit/evidence')
from flint import ctx
import galerkin_fast as gf
for c, N, prec in [(2, 60, 200), (3, 60, 200), (4, 80, 256), (5, 100, 320)]:
    ctx.prec = prec
    M = gf.build(c, N)
    lam, hist, x = gf.lam_min_inv(M, 4)
    print(c, N, lam.mid().str(8, radius=False), flush=True)
