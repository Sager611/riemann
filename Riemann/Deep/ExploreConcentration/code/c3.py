# explore44 c3: values of zeta and DH on [T, 2T] at sigma; Littlewood means, tails vs Bohr lift (c2) and Herbst; Littlewood check with known DH zeros.
import numpy as np, mpmath as mp, json, time, sys
sys.path.insert(0, '/Users/adriansager/playground/riemann/.scratch/explore44')
from dhlib import hurwitz, f as fdh
t0 = time.time()
T = float(sys.argv[1]); h = float(sys.argv[2]); sigs = [float(v) for v in sys.argv[3].split(',')]
mp.mp.dps = 20
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
def fmp(s): return mp.power(5, -s)*(mp.zeta(s, mp.mpf(1)/5) + kap*mp.zeta(s, mp.mpf(2)/5) - kap*mp.zeta(s, mp.mpf(3)/5) - mp.zeta(s, mp.mpf(4)/5))
M = int(2*T/2) + 60
for s in [0.6 + 1j*(T + 0.37), 0.7 + 1j*(2*T - 0.11)]:
    a = hurwitz(np.array([s]), 1.0, M)[0]; b = complex(mp.zeta(mp.mpc(s.real, s.imag)))
    c = fdh(np.array([s]), M=M)[0]; d = complex(fmp(mp.mpc(s.real, s.imag)))
    print('check zeta', s, abs(a - b)/abs(b), ' DH', abs(c - d)/abs(d), flush=True)
t = np.arange(T, 2*T, h)
res = {'T': T, 'h': h}
for sig in sigs:
    s = sig + 1j*t
    z = np.empty(len(t), dtype=complex); F = np.empty(len(t), dtype=complex)
    for i in range(0, len(t), 4096):
        z[i:i+4096] = hurwitz(s[i:i+4096], 1.0, M); F[i:i+4096] = fdh(s[i:i+4096], M=M)
    lz = np.log(np.abs(z)); lf = np.log(np.abs(F))
    rs = [0.5, 1, 1.5, 2, 3, 4, 5, 6]
    r = dict(zeta_mean=float(lz.mean()), zeta_var=float(lz.var()), dh_mean=float(lf.mean()), dh_var=float(lf.var()),
             zeta_low={str(q): float(np.mean(lz < -q)) for q in rs}, dh_low={str(q): float(np.mean(lf < -q)) for q in rs},
             zeta_min=float(lz.min()), dh_min=float(lf.min()), zeta_up={str(q): float(np.mean(lz > q)) for q in rs})
    # local minima of log|DH| below -2: candidate off-line zeros close to the line Re s = sigma
    idx = np.nonzero((lf[1:-1] < lf[:-2]) & (lf[1:-1] < lf[2:]) & (lf[1:-1] < -2.0))[0] + 1
    r['dh_deep_minima'] = [(float(t[i]), float(lf[i])) for i in idx[:40]]
    idz = np.nonzero((lz[1:-1] < lz[:-2]) & (lz[1:-1] < lz[2:]) & (lz[1:-1] < -2.0))[0] + 1
    r['zeta_deep_minima'] = [(float(t[i]), float(lz[i])) for i in idz[:40]]
    res[str(sig)] = r
    np.save('/Users/adriansager/playground/riemann/.scratch/explore44/vals_T%g_s%g.npy' % (T, sig), np.stack([lz, lf]).astype(np.float32))
    print('sigma', sig, 'zeta mean %.5f var %.4f min %.2f | DH mean %.5f var %.4f min %.2f | #deep minima zeta %d DH %d' % (r['zeta_mean'], r['zeta_var'], r['zeta_min'], r['dh_mean'], r['dh_var'], r['dh_min'], len(idz), len(idx)), round(time.time() - t0, 1), flush=True)
    print('   P_T(log|zeta|<-r)', {k: round(v, 5) for k, v in r['zeta_low'].items()})
    print('   P_T(log|DH|<-r)  ', {k: round(v, 5) for k, v in r['dh_low'].items()})
json.dump(res, open('/Users/adriansager/playground/riemann/.scratch/explore44/c3_T%g.json' % T, 'w'))
print('done', round(time.time() - t0, 1))
