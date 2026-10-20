import sys
sys.path.insert(0, '/Users/adriansager/playground/riemann/.scratch/weil')
from flint import arb, arb_mat, ctx
import weil
from refute import coeffs
ctx.prec = 384
for (k, N) in ((7, 80), (7, 120), (6, 120)):
    L = arb(k+1).log(); a = L/2; ak = arb(k).log()/2
    M = weil.weil_matrix(L, N); odd = list(range(1, N, 2)); Mo = weil.sub(M, odd)
    c = coeffs(N, a, ak, a, lambda y: (1 if y > 0 else -1)*(-abs(y)/2).exp())
    nrm2 = 2*(arb(k)**-0.5 - arb(k+1)**-0.5)
    cv = arb_mat(len(odd), 1)
    for i, j in enumerate(odd): cv[i, 0] = c[j]/nrm2.sqrt()
    q = (cv.transpose()*Mo*cv)[0, 0]; m = (cv.transpose()*cv)[0, 0]
    print(f"window log{k+1} N={N}: Q(P_N e_{k}) = {q.mid().str(8, radius=False)}, ||P_N e||^2 = {m.mid().str(6, radius=False)}, lam_min odd = {weil.lam_min(Mo).str(6, radius=False)}", flush=True)
