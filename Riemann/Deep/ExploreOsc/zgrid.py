# zgrid.py (explore60): Z on t in [0, TMAX], step H, for zeta, dh, chi8, chim8; saves logs/grid_<name>.npz; spot checks vs mpmath.
import sys, time, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreOsc')
import oscz as Z
TMAX = float(sys.argv[1]); H = float(sys.argv[2]); names = sys.argv[3].split(',')
t = np.round(np.arange(0.0, TMAX + H/2, H), 10)
for nm in names:
    t0 = time.time()
    z, im = Z.Z(nm, t)
    chk = []
    for tc in (10.5, 100.3, 0.5*TMAX + 0.37, TMAX - 0.61):
        a = Z.Z(nm, np.array([tc]))[0][0]; b = Z.Z_mp(nm, tc); chk.append(abs(a - b))
    sc = int((np.sign(z[:-1])*np.sign(z[1:]) < 0).sum())
    np.savez_compressed('Riemann/Deep/ExploreOsc/logs/grid_%s.npz' % nm, t=t, z=z)
    print('%s: n=%d sign changes=%d max|Im|/max|Z|=%.1e max|Z-Zmp|=%.1e Z(0)=%.6f %.1fs' %
          (nm, len(t), sc, np.abs(im).max()/np.abs(z).max(), max(chk), z[0], time.time() - t0), flush=True)
