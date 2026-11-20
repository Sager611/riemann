# defl_falsify.py (wsgsaudit): falsification of WeilSimple's NEW code path (the rank-one deflation by g, cert_simple.py sector defl).
# Mode 'orth': replace g by a vector orthogonal to the truncated ground vector xi_N. Then f = xi_N has <f, g> = 0 and
# Q(f) = e1(N) <= U < eps, so Q(f) - eps||f||^2 + M<f,g>^2 < 0 and the certificate MUST fail. Mode 'normal': unchanged (control).
# Usage: defl_falsify.py mode c eps bits K2   (env WS_NG, WS_M, WZ_RSCAN, WZ_SCAN as for cert_simple.py)
import sys
sys.path.insert(0, '/tmp/wsgs_copy')
import cert_simple as CS
from flint import arb
mode = sys.argv[1]
orig = CS.ground_vec
def gv(T, n, its=8):
    xi = orig(T, n, its)
    if mode == 'normal': return xi
    e = [arb(1) if i == 0 else arb(0) for i in range(n)]
    d = sum((e[i]*xi[i] for i in range(n)), arb(0)); nn = sum((x*x for x in xi), arb(0))
    g = [(e[i] - d/nn*xi[i]).mid() for i in range(n)]
    print('orth check <g, xi> =', sum((g[i]*xi[i] for i in range(n)), arb(0)).str(5), file=sys.stderr)
    return g
CS.ground_vec = gv
sys.argv = ['cert_simple.py', sys.argv[2], 'defl'] + sys.argv[3:]
CS.main()
