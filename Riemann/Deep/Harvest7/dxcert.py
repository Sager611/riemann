# Arb sign certificate: printed Dimitrov-Xu transform F[Phi_{2,y}](x) = 2 Re(Xi'(z)^2 - Xi(z)Xi''(z)), z = x+iy,
# Xi(z) = xi(1/2+iz) => Xi'^2 - Xi Xi'' = xi(s)xi''(s) - xi'(s)^2, s = 1/2 - y + ix.  Opposite certified signs at two x
# (same y in (0,1/2)) => a real zero (continuity) => translates of Phi_{2,y} not dense in L^1 (Wiener).
# Also prints Jensen's corrected J = |Xi'|^2 - Re(Xi conj Xi'') = |xi'|^2 + Re(xi conj xi'') (>= 0 under RH).
import sys
from flint import acb, arb, acb_series, ctx
ctx.prec = 160
def jet(s0):
    t = acb_series([s0, 1], prec=3)
    f = (t * (t - 1) / 2) * ((-t / 2) * acb.pi().log() + (t / 2).lgamma()).exp() * t.zeta()
    return f[0], f[1], 2 * f[2]
y = arb(sys.argv[1]); 
for xs in sys.argv[2:]:
    x = arb(xs)
    s0 = acb(arb(1) / 2 - y, x)
    a, b, c = jet(s0)
    R = (a * c - b * b).real
    J = (b * b.conjugate()).real + (a * c.conjugate()).real
    n = (a * a.conjugate()).real
    sg = '+' if R > 0 else ('-' if R < 0 else '?')
    print('y=%s x=%s  sign(R)=%s  R/|xi|^2=%s  J/|xi|^2=%s' % (sys.argv[1], xs, sg, (R / n).str(5, radius=False), (J / n).str(5, radius=False)), flush=True)
