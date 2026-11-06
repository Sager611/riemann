import sys, time; sys.path.insert(0, __file__.rsplit('/', 1)[0])
# adversarial (b): xi-like kernel with gamma1, gamma2 merged into gamma_m +- i delta (see offline_np.py)
import numpy as np, mpmath as mp
from offline_np import *
mp.mp.dps = 20
def xi(s): return mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
for delta in [0.05, 0.25, 0.45]:
    u, F0, F1, F2, (gm, d, A, B, N, dp10, dp20) = make(delta)
    print('delta %.2f: gamma_m %.10f A %.6f B %.6f psi1\'(0) %.1e psi2\'(0) %.1e' % (delta, gm, A, B, dp10, dp20))
    h = u[1] - u[0]
    def tr(z):
        f = F0*np.cos(z*u); return 4*h*(f.sum() - 0.5*f[0] - 0.5*f[-1])
    for z in [3.0, 8.0, 12.0, 25.0, complex(gm, d)]:
        zz = mp.mpc(z); XR = complex(xi(mp.mpf(0.5) + 1j*zz)*N(zz**2)/((zz**2 - G1**2)*(zz**2 - G2**2)))
        k = complex(tr(z))
        print('   z=%-28s kernel transform %-45s Xi*R %s' % (z, k, XR))
    sel = (u > 0) & (u <= 2.4)
    with np.errstate(divide='ignore', invalid='ignore'):
        lc = (F0*F2 - F1**2)/F0**2; dl = F1/F0
    print('   K(0) %.10f (log K)\'\'(0) %.6f t_b %.4f | on (0,2.4]: min K %.3e, max K\'/K %.5f at %.3f, max (log K)\'\' %.5f at %.3f' % (
        F0[0], F2[0]/F0[0], -2*F2[0]/F0[0], F0[sel].min(), dl[sel].max(), u[sel][dl[sel].argmax()], lc[sel].max(), u[sel][lc[sel].argmax()]))
    # tail check u in (2.4, 3.2]: relative size of the psi correction
    tail = (u > 2.4) & (u <= 3.2) & (F0 > 0)
    P0 = Phi_np(u[tail])
    print('   tail (2.4,3.2]: max |Phi_off/Phi - 1| = %.2e' % np.max(np.abs(F0[tail]/P0 - 1)))
