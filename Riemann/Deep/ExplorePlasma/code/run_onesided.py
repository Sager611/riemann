import sys, time; sys.path.insert(0, __file__.rsplit('/', 1)[0])
# T1: one-sided transforms Xi+(z) = 2 int_0^inf K e^{izu} du for K = Phi, PhiDH: no zeros on Im z >= 0 (Enestrom-Kakeya), zeros in Im z < 0
import numpy as np, mpmath as mp
from onesided import *
t0 = time.time()
for K, nm, box in [(Phi_np, 'zeta', (100, -30, 401, 240)), (PhiDH_np, 'DH', (80, -14, 321, 112))]:
    F = make_onesided(K)
    print(nm, 'Xi+(0) =', F(0.0)[0].real, ' Xi(5) = Xi+(5)+Xi+(-5) =', (F(5.0) + F(-5.0))[0].real, flush=True)
    xs = np.linspace(-60, 60, 241); ys = np.linspace(0, 20, 41)
    X, Y = np.meshgrid(xs, ys); V = np.abs(F((X + 1j*Y).ravel()))
    print(nm, 'min |Xi+| on the grid [-60,60]x[0,20] (upper half-plane):', V.min(), flush=True)
    zs = find_zeros(F, np.linspace(0.0, box[0], box[2]), np.linspace(box[1], -0.05, box[3]))
    print(nm, 'zeros of Xi+ with 0<Re<%d, %d<Im<0: %d' % (box[0], box[1], len(zs)))
    for z in zs: print('   %.8f %+.8fi' % (z.real, z.imag))
    print('   S(x) = Im Xi+(x) > 0 on (0,60]?', bool(np.all(F(np.linspace(0.01, 60, 6000)).imag > 0)), flush=True)
    if nm == 'zeta':
        from kern import ev, TZ
        mp.mp.dps = 30
        g = lambda z: 2*mp.quad(lambda u: ev(TZ, u)*mp.e**(1j*z*u), mp.linspace(0, 4.5, 19), method='gauss-legendre')
        r = mp.findroot(g, mp.mpc(zs[0].real, zs[0].imag))
        print('   mpmath (dps 30) refine of the first zero:', mp.nstr(r, 15), ' |Xi+| =', mp.nstr(abs(g(r)), 3), flush=True)
print('time %.1f' % (time.time() - t0))
