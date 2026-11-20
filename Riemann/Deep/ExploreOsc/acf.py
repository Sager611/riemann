# acf.py (explore60): Z in the class S? Normalised unfolded autocorrelation (1/V) int F(v)F(v+u) dv vs sinc(pi u), and the
# time-domain Rice ratio (1/pi) sqrt(int Z'^2 / int Z^2) / (theta'/pi) (flat-measure value 1/sqrt 3) vs the actual density.
import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreOsc')
import oscz as O
h = 0.01
for nm in ('zeta', 'dh', 'chi8', 'chim8'):
    g = np.load('Riemann/Deep/ExploreOsc/logs/grid_%s.npz' % nm); t, z = g['t'], g['z']
    m = (t >= 1000) & (t <= 2000)
    tt, zz = t[m], z[m]
    th = O.theta(nm, tt); om = np.gradient(th, h)
    zp = np.gradient(zz, h)
    sc = int((np.sign(zz[:-1])*np.sign(zz[1:]) < 0).sum())
    gc = (th[-1] - th[0])/np.pi
    rice_loc = np.sqrt(np.sum(zp**2/om**2)/np.sum(zz**2))      # local units: divide Z' by theta'
    var = np.mean(zz**2/(2*om))
    # unfolded grid v = theta/pi, uniform step dv; F = Z/sqrt(2 theta')
    v = th/np.pi; dv = 0.005; vg = np.arange(v[0], v[-1], dv)
    F = np.interp(vg, v, zz/np.sqrt(2*om))
    row = []
    for u in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
        k = int(round(u/dv)); c = np.mean(F[:-k]*F[k:])/np.mean(F*F)
        row.append('%.2f:%+.3f(%+.3f)' % (u, c, np.sinc(u)))
    print('%-5s [1000,2000] sign changes %d / theta-count %.1f = %.4f | Rice (time domain) %.4f (flat 0.5774) | '
          'E[Z^2/(2theta\')] %.3f | acf u:emp(sinc) %s' % (nm, sc, gc, sc/gc, rice_loc, var, ' '.join(row)), flush=True)
