# ceilaudit7 B: zeta/CUE shift-dependent fourth-order mean G(x) = <Z(0)^2 Z(x)^2>/<Z^4> (normalised units) vs Gaussian block mixtures (G >= 1/3)
import numpy as np, sys
sys.argv = ['x', 'none']
exec(open('b4.py').read().split("mode = sys.argv[1]")[0])
for N in (400, 800):
    l0, _ = F(N, 0.0)
    out = []
    for x in (0.5, 1, 2, 4, 8, 16):
        lx, s = F(N, 2*np.pi*x/N); out.append(f"G({x})={s.real*np.exp(lx-l0):.4f}")
    print(f"N={N} " + " ".join(out) + f"; 1/(2x^2) at x=4,8,16: {1/32:.4f} {1/128:.4f} {1/512:.4f}; 1/M4_CUE={12*(N+1)/((N+2)**2*(N+3)):.2e}")
