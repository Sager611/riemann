# C1: a-point (Levinson-Montgomery-type) identity, smallest case t in [T1,T2].
# V_c = F' - (1/2+c)(chi'/chi) F.  Claim: #zeros(V_c, sigma>1/2) = N+ (0<c) and N+ + N0bar (-1/2<c<0), up to O(1) edge terms.
import numpy as np, math, sys
from fz import Fn, winding, rect
T1, T2 = float(sys.argv[1]), float(sys.argv[2]); kinds = sys.argv[3].split(',')
cs = [3, 1, 0.5, 0.25, 0.1, 0.05, -0.05, -0.1, -0.25, -0.4, -0.49]
for kind in kinds:
    f = Fn(kind, T2 + 10)
    t = np.arange(T1, T2, 0.002); F, _ = f(0.5 + 1j*t); Z = np.real(np.exp(1j*f.theta(t))*F)
    n0 = int(np.sum(np.sign(Z[1:]) != np.sign(Z[:-1])))
    P = rect(-1.5, 2.5, T1, T2, 0.004, 0.02, 0.002); Fb, _ = f(P); ntot, st = winding(Fb, P)
    npl = (ntot - n0)/2
    out = [f"{kind}: [{T1},{T2}] line zeros (Z sign changes) N0bar={n0}; all zeros (winding on [-1.5,2.5]) {ntot:.3f} (max step {st:.2f}); N+ = {npl:.2f}"]
    for eps in (0.005, 0.02):
        P = rect(0.5+eps, 3, T1, T2, 0.0005, 0.02, 0.0005); Fb, _ = f(P); w, st = winding(Fb, P)
        out.append(f"  zeros of F with sigma>1/2+{eps}: {w:.3f} (max step {st:.2f})")
    P = rect(0.5, 3, T1, T2, 0.002, 0.02, 0.001); Fb, Fpb = f(P); ch = f.chilog(P)
    row = []
    for c in cs:
        w, st = winding(Fpb - (0.5 + c)*ch*Fb, P); row.append(f"c={c:+}: {w:.2f}" + ("" if st < 1.5 else f"(step {st:.1f}!)"))
    out.append("  #zeros of V_c right of 1/2: " + "; ".join(row))
    print("\n".join(out), flush=True)
