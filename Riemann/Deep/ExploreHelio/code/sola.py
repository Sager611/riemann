# sola.py (explore23, T3): SOLA/OLA kernels as Weil test functions. Q(f) = sum_rho F(g_rho) F(-g_rho), F(z) = int_W f(x) e^{izx} dx.
# An off-line pair (gbar, delta) contributes 2 Re A(gbar - i delta) = 2A(gbar) - delta^2 A''(gbar) + O(delta^4), A = F F*.
import json, numpy as np
PI = np.pi
gz = np.load('zeta_zeros.npy')
GA, GB = gz[5], gz[6]; GBAR = 0.5*(GA + GB)
def FT(fx, xs, z):
    # F(z) = int f(x) e^{izx} dx on grid xs (vector z, complex allowed), trapezoid
    dx = xs[1] - xs[0]; wts = np.full(len(xs), dx); wts[0] = wts[-1] = dx/2
    out = np.empty(len(z), dtype=complex)
    for i in range(0, len(z), 2000):
        out[i:i+2000] = (np.exp(1j*np.outer(z[i:i+2000], xs))*(fx*wts)[None, :]).sum(1)
    return out
res = []
for L in [1.0, 2.0, 3.0, 4.0]:
    xs = np.linspace(-L/2, L/2, 801)
    for kind in ['sola', 'dipole']:
        for t0 in [GBAR, GBAR + 1.0, GBAR + 3.0]:
            if kind == 'sola':
                fx = np.cos(PI*xs/L)*np.cos(t0*xs)            # Backus-Gilbert spread-optimal taper, peak at +-t0
            else:
                fx = np.sin(2*PI*xs/L)*np.sin(t0*xs)          # odd taper x odd carrier: simple node at +-t0
            for delta in [0.15, 0.3]:
                zon = np.concatenate([gz[:5], gz[7:]])            # zeta zeros minus 6,7
                Fon = FT(fx, xs, zon).real
                Qbase = 2*np.sum(Fon**2)                         # +-gamma
                zc = np.array([GBAR - 1j*delta, GBAR + 1j*delta])
                Fp = FT(fx, xs, zc); Fm = FT(fx, xs, -zc)
                off = 2*np.sum(Fp*Fm).real                       # both members, +- (factor 2 for -gamma side)
                Ford = FT(fx, xs, np.array([GBAR])).real[0]
                dbl = 4*Ford**2                                   # double ordinate at gbar, +- sides
                Fz = FT(fx, xs, gz[5:7]).real; Qzeta = Qbase + 2*np.sum(Fz**2)
                # curvature: A(t) = F(t)^2 (f real even -> F real even), A'' by finite differences
                h = 1e-3; Ft = FT(fx, xs, np.array([GBAR - h, GBAR, GBAR + h])).real; A = Ft**2
                App = (A[0] - 2*A[1] + A[2])/h**2
                res.append({'L': L, 'kind': kind, 't0-gbar': round(t0 - GBAR, 3), 'delta': delta, 'Q_zeta': float(Qzeta), 'Q_mv': float(Qbase + off),
                            'Q_mv_ordinates(double)': float(Qbase + dbl), 'offline_minus_double': float(off - dbl), 'pred_-2delta^2A\'\'': float(-2*delta**2*App)})
json.dump(res, open('sola.json', 'w'), indent=1)
for r in res:
    print(r['L'], r['kind'], r['t0-gbar'], r['delta'], '%.4g %.4g %.4g %.4g %.4g' % (r['Q_zeta'], r['Q_mv'], r['Q_mv_ordinates(double)'], r['offline_minus_double'], r["pred_-2delta^2A''"]))
