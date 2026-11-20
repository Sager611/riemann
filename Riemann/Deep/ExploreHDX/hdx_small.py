# hdx_small.py (explore64): short links. e1(c'), e2(c') (float Galerkin, Legendre, both parities) against the Kronecker ceiling
# Lam(c, c') (hdx_ceiling.optimize, float) for c = 10 and 47: where can a scalar (Garland) or rank-1 (Oppenheim) link floor pass the ceiling?
import sys, math, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreFunctional/code'); sys.path.insert(0, 'Riemann/Deep/ExploreHDX')
from xform import Form
from hdx_ceiling import optimize
from flint import fmpq

def e12(c, n=24):
    ev = []
    for sec in ['even', 'odd']:
        Q, _ = Form(c).Q(n, sec); ev += list(np.linalg.eigvalsh((Q + Q.T)/2)[:3])
    ev.sort(); return ev[0], ev[1]

if __name__ == '__main__':
    print("# c_link, e1, e2 (float Galerkin n=24/sector), Lam(10, c_link), Lam(47, c_link) (float ceiling, N=1500)")
    for cl in [1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0, 2.5, 3.162]:
        e1, e2 = e12(cl)
        clq = fmpq(int(round(cl*10**6)), 10**6)
        L10 = optimize(fmpq(10), clq, 1500)[1]; L47 = optimize(fmpq(47), clq, 1500)[1]
        print(f'{cl:7.3f}  e1 {e1: .4f}  e2 {e2: .4f}  Lam10 {L10:.4f}  Lam47 {L47:.4f}', flush=True)
