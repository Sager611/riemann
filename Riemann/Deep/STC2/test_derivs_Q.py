# test_derivs_Q.py -- check certK2.py tail form Qv against t^2 U - 4d^3 u - K5|u|^3 and Qder against finite differences (D6000)
import sys
sys.argv = ['x', 'D6000', 'box']
exec(open('Riemann/Deep/STC2/certK2.py').read().split("t0 = time.time()")[0])
from flint import arb, acb
for (b, u, sgn) in ((0.3, 0.007, 1), (0.45, -0.004, -1), (0.1, 0.009, 1)):
    t = 1/u; s = 1 if u > 0 else -1
    q = Qv(arb(b), arb(u), sgn)
    # t^2 U(b,t) directly from vg (s = +1 coefficients, t may be negative): use H + f
    U = arb(0)
    for x in (d + arb(b), d - arb(b)):
        z = acb(x, arb(t))
        for (X, j), c in zip(B, C):
            w = z**(-j); U += c*(w.real if X == 'R' else -w.imag)
    q2 = arb(t)**2*U - 4*d**3*arb(u) - K5*sgn*arb(u)**3
    e = 1e-5
    g1, g2 = Qder(arb(b), arb(u))
    n1 = (Qv(arb(b+e), arb(u), sgn) - Qv(arb(b-e), arb(u), sgn))/(2*e)
    n2 = (Qv(arb(b+e), arb(u), sgn) - 2*Qv(arb(b), arb(u), sgn) + Qv(arb(b-e), arb(u), sgn))/(e*e)
    print(b, u, 'Q %.12e vs %.12e' % (float(q.mid()), float(q2.mid())), 'Qb %.8e/%.8e  Qbb %.6e/%.6e' % (float(g1.mid()), float(n1.mid()), float(g2.mid()), float(n2.mid())))

