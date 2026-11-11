
# T3: "black hole in a box": mirror at height a in the cusp -> normal modes = zeros of a^s X(2s) + a^{1-s} X(2s-1).
# Hermite-Biehler: all on Re s = 1/2 if g_X(u) = Re(xi_X'/xi_X)(1+iu) > 0 for all u (a >= 1), i.e. Lagarias positivity at sigma = 1.
import numpy as np
from flint import acb, arb, ctx, dirichlet_char
ctx.prec = 128  # derivative by central differences, h = 1e-15
PI = acb(arb.pi())
h = acb(1e-15)
def xi_zeta(w):
    return w*(w-1)/2 * (PI**(-w/2)) * (w/2).gamma() * w.zeta()
g1, g2 = 14.134725141734693790, 21.022039638771554993
gb, dl = (g1+g2)/2, 0.3
def Rmod(w):
    u = w - acb(0.5)
    num = ((u - acb(0, gb))**2 - dl*dl) * ((u + acb(0, gb))**2 - dl*dl)
    den = (u*u + g1*g1) * (u*u + g2*g2)
    return num/den
def xi_mod(w): return xi_zeta(w) * Rmod(w)
ca, cb = dirichlet_char(5, 2), dirichlet_char(5, 3)
def gam5(s): return (acb(5)/PI)**(s/2) * ((s+1)/2).gamma()
t0 = 10.0
s0 = acb(0.5, t0)
e0 = gam5(s0); e0 = e0/abs(e0)
ka = acb(ca.hardy_z(acb(t0)).real) / (e0 * ca.l(s0))
kb = acb(cb.hardy_z(acb(t0)).real) / (e0 * cb.l(s0))
def fDH(w): return ka*ca.l(w) + kb*cb.l(w)
def xi_DH(w): return gam5(w) * fDH(w)
# DH zero calibration: see t3_dhcal.py
zz = (gam5(acb(0.5,50.0))/abs(gam5(acb(0.5,50.0)))) * fDH(acb(0.5,50.0))
print("reality of e^{i theta} f_DH on the line at t=50: Im =", zz.imag.mid().str(3), " Re =", zz.real.mid().str(8), " ZA+ZB =", (ca.hardy_z(acb(50.0)).real + cb.hardy_z(acb(50.0)).real).mid().str(8))
def logder(F, w): return (F(w + h) - F(w - h)) / (2*h) / F(w)
res = {}
for name, F, U in (("zeta", xi_zeta, 200.0), ("moved pair", xi_mod, 200.0), ("DH", xi_DH, 200.0)):
    us = np.arange(0.05, U, 0.05)
    g = np.array([float(logder(F, acb(1.0, float(u))).real.mid()) for u in us])
    k = int(np.argmin(g))
    res[name] = (us, g)
    print("%-10s min_u Re(xi'/xi)(1+iu) over (0,%g]: %.6f at u = %.2f ; mean %.4f" % (name, U, g[k], us[k], g.mean()))
# box spectrum: b(t) = Re(a^{it} xi(1+2it)), zeros vs phase count, t in (0, 100]
ts = np.arange(0.02, 100.0, 0.02)
for name, F in (("zeta", xi_zeta), ("moved pair", xi_mod), ("DH", xi_DH)):
    X = []
    for t in ts:
        v = F(acb(1.0, 2*float(t))); X.append(complex(float(v.real.mid()), float(v.imag.mid())))
    X = np.array(X); argX = np.unwrap(np.angle(X))
    for a in (1.0, 7.0555, 20.0):
        ph = argX + ts*np.log(a); b = np.abs(X)*np.cos(ph)
        sc = int(np.sum(np.abs(np.diff(np.signbit(b).astype(int)))))
        print("box %-10s a=%-7g sign changes of Re E on (0,100]: %d ; phase increase/pi: %.3f ; min dphase/dt %.4f" % (name, a, sc, (ph[-1]-ph[0])/np.pi, np.min(np.diff(ph))/0.02))
