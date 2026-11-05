import sys, math
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ChekAudit/evidence')
from flint import ctx
import galerkin_fast as gf
for c, N, prec in [(2, 20, 200), (2, 40, 200), (2, 80, 200), (2, 160, 256), (3, 120, 256)]:
    ctx.prec = prec
    M = gf.build(c, N)
    lam, hist, x = gf.lam_min_inv(M, 5)
    print(c, N, lam.mid().str(8, radius=False), [h[1].str(2, radius=False) for h in hist[-2:]], flush=True)
