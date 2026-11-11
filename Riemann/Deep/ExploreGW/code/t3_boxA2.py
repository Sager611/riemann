
import numpy as np
from flint import acb, arb, ctx, dirichlet_char
ctx.prec = 96
PI = acb(arb.pi())
ch = {n: dirichlet_char(13, n) for n in (3, 9)}
def gam(s): return (acb(13)/PI)**(s/2) * (s/2).gamma()
def mid(z): return acb(z.real.mid(), z.imag.mid())
def kappa(n, t0=10.0):
    s0 = acb(0.5, t0); e0 = gam(s0); e0 = e0/abs(e0)
    return mid(acb(ch[n].hardy_z(acb(t0)).real) / (e0 * ch[n].l(s0)))
K = {n: kappa(n) for n in ch}
def xiA(w): return gam(w) * (K[3]*ch[3].l(w) + K[9]*ch[9].l(w))
def conj(z): return acb(z.real, -z.imag)
def B(s):  # E(s) + E#(s), E(s) = xi_A(2s), E#(s) = conj(E(1 - conj s))
    return xiA(2*s) + conj(xiA(2*(1 - conj(s))))
h = acb(1e-12)
found = []
for z0 in [(0.52, 38.55), (0.55, 38.60), (0.60, 38.60), (0.65, 38.65), (0.58, 38.5), (0.53, 38.7)]:
    s = acb(*z0)
    for _ in range(60):
        fs = mid(B(s)); d = mid((mid(B(s+h)) - mid(B(s-h)))/(2*h)); st = fs/d
        if float(abs(st).mid()) > 0.03: st = st*(0.03/float(abs(st).mid()))
        s = mid(s - st)
        if float(abs(st).mid()) < 1e-14: break
    found.append((round(float(s.real.mid()), 8), round(float(s.imag.mid()), 6), "%.1e" % float(abs(B(s)).mid())))
print("Newton zeros of the block-A box function near rho/2 = 0.5915+38.596i:", found)
