# pyright: reportMissingImports=false
"""Desogus Thm 8.6 test: T_k >= g |e_k><e_k| (g >= 1.878) <=> Q(v) >= g |<e_k,v>|^2 for all odd v on (-a_{k+1},a_{k+1}).
inf_v Q(v)/|<e,v>|^2 = 1/<e, A_{k+1}^{-1} e>. Galerkin (odd Legendre) gives an UPPER bound on that infimum."""
import sys, time
sys.path.insert(0, '/Users/adriansager/playground/riemann/.scratch/weil')
import flint
from flint import arb, arb_mat, ctx
import weil

def coeffs(N, a, lo, hi, prof, nq=80):
    """c_j = int_{lo<|y|<hi} prof(y) phi_j(y) dy, phi_j = leg_j(y/a)/sqrt(a); prof odd."""
    c = [arb(0)]*N
    half = (hi - lo)/2; mid = (hi + lo)/2
    for (x, w) in weil.gauss(nq):
        for s in (1, -1):
            y = s*(mid + half*x)
            u = weil.leg_on(N, y/a)
            p = prof(y)*w*half/a.sqrt()
            c = [c[j] + u[j]*p for j in range(N)]
    return c

def run(k, N, prec):
    ctx.prec = prec
    t0 = time.time()
    L = arb(k+1).log(); a = L/2
    ak = arb(k).log()/2
    M = weil.weil_matrix(L, N)
    odd = list(range(1, N, 2))
    Mo = weil.sub(M, odd)
    profs = {
      'ground sgn(y)e^{-|y|/2} (paper e_k)': lambda y: (1 if y > 0 else -1)*(-abs(y)/2).exp(),
      'polar t_k = sinh(y/2) on collar': lambda y: (y/2).sinh(),
      'sign indicator of collar': lambda y: arb(1 if y > 0 else -1),
    }
    lam = weil.lam_min(Mo)
    print(f"k={k}: window log{k+1}, collar ({float(ak):.5f},{float(a):.5f}), N={N}, prec={prec}, lam_min(A_{k+1},odd Galerkin)={lam.str(6, radius=False)}", flush=True)
    for name, pf in profs.items():
        c = coeffs(N, a, ak, a, pf)
        nrm2 = sum((c[j]*c[j] for j in odd), arb(0))
        # exact L2 norm of profile on collar via quadrature
        tot = arb(0); half = (a-ak)/2; mid = (a+ak)/2
        for (x, w) in weil.gauss(80):
            y = mid + half*x; tot += 2*pf(y)**2*w*half
        cv = arb_mat(len(odd), 1)
        for i, j in enumerate(odd): cv[i, 0] = c[j]/tot.sqrt()
        x = Mo.solve(cv)
        q = (cv.transpose()*x)[0, 0]
        print(f"  e = {name}: captured L2 mass {float(nrm2/tot):.6f}; <e,A^-1 e> >= {q.mid().str(6, radius=False)}; inf Q(v)/|<e,v>|^2 <= {(1/q).mid().str(6, radius=False)}  (paper claims >= g_k >= 1.878)", flush=True)
        if name.startswith('ground'):
            # explicit trial v = A^-1 e (Galerkin), report Q(v), <e,v>
            Qv = (x.transpose()*Mo*x)[0, 0]; ev = (cv.transpose()*x)[0, 0]
            print(f"    trial v=M^-1 e: Q(v)={Qv.mid().str(6, radius=False)} <e,v>={ev.mid().str(6, radius=False)}", flush=True)
    print(f"  t={time.time()-t0:.1f}s", flush=True)

if __name__ == '__main__':
    N, prec = int(sys.argv[1]), int(sys.argv[2])
    for k in sys.argv[3:]: run(int(k), N, prec)
