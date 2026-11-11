# Certify that Z(t) has no zero on [g_{n-1}, g_{n+1}] (zero-free Gram block): adaptive bisection with Arb enclosures of Z over
# subintervals (t as a ball); a piece is accepted when its enclosure excludes 0.  usage: certify_zerofree.py n
import sys
from flint import arb, acb, ctx
import mpmath
ctx.prec = 128
PI = arb.pi()
def theta(t): return acb(arb(0.25), t/2).lgamma().imag - t/2*PI.log()
def Zball(a, b):
    t = arb((a + b)/2, (b - a)/2)
    th = theta(t)
    return (acb(th.cos(), th.sin())*acb(arb(0.5), t).zeta()).real
n = int(sys.argv[1]); mpmath.mp.dps = 30
lo = mpmath.grampoint(n - 1); hi = mpmath.grampoint(n + 1)
a0 = arb(mpmath.nstr(lo, 25)) - arb('1e-15'); b0 = arb(mpmath.nstr(hi, 25)) + arb('1e-15')   # slightly wider than the block
stack = [(a0.mid(), b0.mid())]; pieces = 0; minabs = None; sign = None
while stack:
    a, b = stack.pop()
    z = Zball(arb(a), arb(b))
    if z > 0 or z < 0:
        s = 1 if z > 0 else -1
        assert sign is None or s == sign, 'sign change: a zero lies in the block'
        sign = s; pieces += 1
        lb = float(abs(z).lower()); minabs = lb if minabs is None else min(minabs, lb)
    else:
        m = (arb(a) + arb(b))/2
        if float((arb(b) - arb(a)).mid()) < 1e-12: raise SystemExit('cannot separate from 0 at ' + str(m))
        stack.append((a, m.mid())); stack.append((m.mid(), b))
print(f'n={n}: Z has constant sign {sign:+d} on [g_(n-1), g_(n+1)] = [{mpmath.nstr(lo, 20)}, {mpmath.nstr(hi, 20)}]; {pieces} certified pieces, min |Z| lower bound {minabs:.3e}: zero-free Gram block (certified)')
