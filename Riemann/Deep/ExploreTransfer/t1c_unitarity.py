# t1c_unitarity.py -- explore3: are eigenvalue branches unimodular along Re b = 1/2 (cusp-form branches) and along
# Re b = 1/4 (zeta-zero branch)?  Track the eigenvalue through +-1 at a cusp form / through 1 at rho_1/2 along the
# vertical line, print |lambda| and arg.
import sys, time
import numpy as np
sys.path.insert(0, '.')
from t1_flint import build_M
from t1b_branches import track

if __name__ == '__main__':
    t0 = time.time()
    cases = [('even cusp r=13.7798, sigma=1/2', 0.5, 13.779751351890738, 1.0),
             ('odd cusp r=9.5337, sigma=1/2', 0.5, 9.533695261353557, -1.0),
             ('zeta rho1/2 t=7.0674, sigma=1/4', 0.25, 14.134725141734693790/2, 1.0),
             ('zeta rho2/2 t=10.5110, sigma=1/4', 0.25, 21.022039638771554993/2, 1.0)]
    for name, s, t, tgt in cases:
        ts = [t - 1.0 + 0.05*i for i in range(41)]
        pts = [complex(s, x) for x in ts]
        lam, ev = track(pts, 100, 20, tgt)
        print(name)
        print('   dt   |lambda|   arg/pi   spectral-radius')
        for i in range(0, 41, 4):
            print('  %+.2f  %.5f  %+.4f   %.4f' % (ts[i] - t, abs(lam[i]), np.angle(lam[i]) / np.pi, abs(ev[i]).max()))
        mods = np.array([abs(l) for l in lam])
        print('  max | |lambda|-1 | over dt in [-1,1]: %.4f' % abs(mods - 1).max())
    print('time %.1fs' % (time.time() - t0))
