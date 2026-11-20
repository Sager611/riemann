# analyze.py (explore60): E1 local power spectra (low-band fraction, Rice ratio), E2 frequency-0 coefficient (int Z),
# E3 sign pattern / mollified Z low band, E4 Kedem higher-order crossings. Reads logs/grid_<name>.npz (zgrid.py).
import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreOsc')
import oscz as O
NAMES = ['zeta', 'dh', 'chi8', 'chim8']
G = {nm: np.load('Riemann/Deep/ExploreOsc/logs/grid_%s.npz' % nm) for nm in NAMES}
h = 0.01


def Omega(nm, t):  # local top frequency theta'(t)
    return (O.theta(nm, t + 1e-4) - O.theta(nm, t - 1e-4))/2e-4


def spec(t, x, tc, L, pad=8):
    m = (t >= tc - L/2) & (t <= tc + L/2)
    xs = x[m]*np.hanning(m.sum())
    n = pad*len(xs)
    P = np.abs(np.fft.rfft(xs, n))**2
    w = 2*np.pi*np.fft.rfftfreq(n, h)
    return w, P


def frac(w, P, cut):
    return P[w < cut].sum()/P.sum()


def E1():
    print('E1 local power spectrum: fraction of power at |w| < a*Omega (flat prediction a), Rice density/(Omega/pi), actual/(Omega/pi)')
    for nm in NAMES:
        t, z = G[nm]['t'], G[nm]['z']
        for tc, L in ((150, 60), (500, 60), (1000, 80), (1800, 80)):
            Om = Omega(nm, tc)
            w, P = spec(t, z, tc, L)
            fr = [frac(w, P, a*Om) for a in (0.1, 0.25, 0.5, 0.75, 1.0, 1.25)]
            rice = np.sqrt((w**2*P).sum()/P.sum())/np.pi
            m = (t >= tc - L/4) & (t <= tc + L/4)
            zs = z[m]
            act = (np.sign(zs[:-1])*np.sign(zs[1:]) < 0).sum()/(L/2)
            print('  %-5s tc=%4d Om=%.3f  F(.1,.25,.5,.75,1,1.25)=%s  rice=%.3f  actual=%.3f' %
                  (nm, tc, Om, ' '.join('%.3f' % f for f in fr), rice/(Om/np.pi), act/(Om/np.pi)))


def E2():
    C8 = 2**1.5*np.pi**0.25/3
    print('E2 frequency-0 coefficient I(T) = int_0^T Z; chi8 theory C T^{3/4}, C = %.6f' % C8)
    for nm in NAMES:
        t, z = G[nm]['t'], G[nm]['z']
        I = np.concatenate([[0.0], np.cumsum((z[1:] + z[:-1])*h/2)])
        row = []
        for T in (250, 500, 1000, 1500, 2000):
            j = int(round(T/h))
            row.append('T=%d I=%.2f I/T^.75=%.4f' % (T, I[j], I[j]/T**0.75))
        sel = t >= 100
        r14 = np.abs(I[sel])/t[sel]**0.25
        r8 = np.abs(I[sel] - (C8*t[sel]**0.75 if nm == 'chi8' else 0))/(t[sel]**0.25*np.log(t[sel]))
        print('  %-5s %s | max_{T in[100,2000]} |I|/T^.25=%.3f  |I - main|/(T^.25 log T) max=%.3f' %
              (nm, '; '.join(row), r14.max(), r8.max()))
    # local means over stationary zones: theta'(t) = log n  <=>  t_n = 2 pi n^2 / q
    for nm in ('zeta', 'chi8', 'chim8', 'dh'):
        q = O.FUN[nm][0]
        t, z = G[nm]['t'], G[nm]['z']
        I = np.concatenate([[0.0], np.cumsum((z[1:] + z[:-1])*h/2)])
        out = []
        for n in range(2, 40):
            tn = 2*np.pi*n*n/q
            if tn + 2*np.sqrt(tn) > 1999 or tn - 2*np.sqrt(tn) < 20: continue
            a, b = int(round((tn - np.sqrt(tn))/h)), int(round((tn + np.sqrt(tn))/h))
            out.append('%d:%+.3f' % (n, (I[b] - I[a])/(t[b] - t[a])*tn**0.25))
        print('  %-5s stationary-zone means x t_n^{1/4} (n:value): %s' % (nm, ' '.join(out[:14])))


def mollifier(t, y):
    mu = {1: 1, 2: -1, 3: -1, 5: -1, 6: 1, 7: -1, 10: 1, 11: -1, 13: -1, 14: 1, 15: 1, 17: -1, 19: -1, 21: 1, 22: 1, 23: -1,
          26: 1, 29: -1}
    M = np.zeros(len(t), dtype=complex)
    for m, u in mu.items():
        if m <= y:
            M += u*(1 - np.log(m)/np.log(y))*m**(-0.5 - 1j*t)
    return M


def E3():
    print('E3 low-band fraction F(a) for x = Z, sign(Z), Z|M|^2 (y = 30), window L = 200')
    for nm in NAMES:
        t, z = G[nm]['t'], G[nm]['z']
        for tc in (1000, 1800):
            Om = Omega(nm, tc)
            res = []
            for lab, x in (('Z', z), ('sgnZ', np.sign(z)), ('Z|M|^2', z*np.abs(mollifier(t, 30))**2)):
                w, P = spec(t, x, tc, 200, pad=2)
                res.append('%s:%s' % (lab, '/'.join('%.3f' % frac(w, P, a*Om) for a in (0.25, 0.5, 0.75))))
            print('  %-5s tc=%4d  %s' % (nm, tc, '  '.join(res)))


def E4():
    print('E4 Kedem higher-order crossings: D_k = sign changes of the k-th difference on [20, 2000] (and on [20, 236.6])')
    for nm in NAMES:
        t, z = G[nm]['t'], G[nm]['z']
        th = O.theta(nm, np.array([20.0, 236.6, 2000.0]))
        for (lo, hi), gp in (((20, 2000), (th[2] - th[0])/np.pi), ((20, 236.6), (th[1] - th[0])/np.pi)):
            m = (t >= lo) & (t <= hi)
            x = z[m]
            D = []
            for k in range(0, 7):
                D.append(int((np.sign(x[:-1])*np.sign(x[1:]) < 0).sum()))
                x = np.diff(x)
            print('  %-5s [%d,%g]: D_0..D_6 = %s ; (theta(hi)-theta(lo))/pi = %.1f' % (nm, lo, hi, D, gp))


if __name__ == '__main__':
    for e in sys.argv[1].split(','):
        globals()[e]()
