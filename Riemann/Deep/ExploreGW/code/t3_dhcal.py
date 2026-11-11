
from flint import acb, arb, ctx, dirichlet_char
ctx.prec = 96
PI = acb(arb.pi())
def mid(z): return acb(z.real.mid(), z.imag.mid())
ca, cb = dirichlet_char(5, 2), dirichlet_char(5, 3)
def gam5(s): return (acb(5)/PI)**(s/2) * ((s+1)/2).gamma()
s0 = acb(0.5, 10.0); e0 = gam5(s0); e0 = e0/abs(e0)
ka = mid(acb(ca.hardy_z(acb(10.0)).real)/(e0*ca.l(s0))); kb = mid(acb(cb.hardy_z(acb(10.0)).real)/(e0*cb.l(s0)))
print("kappa_a, kappa_b:", ka.str(8), kb.str(8))
def f(w): return ka*ca.l(w) + kb*cb.l(w)
h = acb(1e-12)
for w0 in ((0.808517, 85.699348), (0.650830, 114.163343), (0.574356, 166.479306), (0.724258, 176.702461)):
    w = acb(*w0)
    for _ in range(30):
        fw = mid(f(w)); d = mid((mid(f(w+h)) - mid(f(w-h)))/(2*h)); st = fw/d; w = mid(w - st)
        if float(abs(st).mid()) < 1e-15: break
    print("DH zero from %s -> %s + %s i, |f| = %s" % (w0, w.real.mid().str(10), w.imag.mid().str(10), abs(f(w)).mid().str(3)))
# classical DH normalisation check: f_class = ((1 - i k)/2) L(chi) + ((1 + i k)/2) L(chibar), k = (sqrt(10-2 sqrt5) - 2)/(sqrt5 - 1)
import math
kk = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
print("chi(2) for label 2:", ca.chi_exponent(2), "of order", ca.order())
w = acb(0.808517, 85.699348)
for lab, A, B in (("chi=label2", ca, cb), ("chi=label3", cb, ca)):
    fc = acb(0.5, -kk/2)*A.l(w) + acb(0.5, kk/2)*B.l(w)
    print("classical DH with", lab, "at 0.808517+85.699348i: |f| =", abs(fc).mid().str(4))
