# Riemann-Siegel Z'(t) at zeros (vectorized) + validation against mpmath. explore55.
import numpy as np, sys
TWOPI = 2*np.pi
def theta(t):
    return t/2*np.log(t/TWOPI) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
def dtheta(t):
    return 0.5*np.log(t/TWOPI) - 1/(48*t**2) - 7/(1920*t**4)
def _c0raw(p):
    return np.cos(2*np.pi*(p*p - p - 1/16))/np.cos(2*np.pi*p)
def c0(p):
    p = np.asarray(p, float); out = _c0raw(p)
    bad = np.abs(np.cos(2*np.pi*p)) < 1e-3
    if bad.any():
        h = 4e-3; out[bad] = 0.5*(_c0raw(p[bad]+h) + _c0raw(p[bad]-h))
    return out
def dc0(p, h=1e-4):
    return (c0(p+h) - c0(p-h))/(2*h)
def zprime(t, chunk=500):
    """Z'(t): derivative of the RS main sum plus derivative of the C0 correction; error O(t^-5/4)-ish for the main part."""
    t = np.asarray(t, float); out = np.empty_like(t)
    for a in range(0, len(t), chunk):
        tt = t[a:a+chunk]; N = np.floor(np.sqrt(tt/TWOPI)).astype(int); Nm = N.max()
        n = np.arange(1, Nm+1, dtype=float); ln = np.log(n)
        th = theta(tt); dth = dtheta(tt)
        ph = th[:, None] - tt[:, None]*ln[None, :]
        w = n**-0.5; mask = (n[None, :] <= N[:, None])
        main = -2*np.sum(mask*w*(dth[:, None] - ln[None, :])*np.sin(ph), axis=1)
        p = np.sqrt(tt/TWOPI) - N; sg = np.where((N-1) % 2 == 0, 1.0, -1.0)
        corr = sg*(-0.25*(tt/TWOPI)**-1.25/TWOPI*c0(p) + (tt/TWOPI)**-0.25*dc0(p)/(2*np.sqrt(TWOPI*tt)))
        out[a:a+chunk] = main + corr
    return out
if __name__ == '__main__':
    import mpmath as mp
    mp.mp.dps = 20
    g = np.loadtxt(sys.argv[1])
    for k in [0, 1, 9, 99, 999, 9999, 49999, 99999]:
        t = g[k]; zp = zprime(np.array([t]))[0]
        ex = float(mp.siegelz(t, derivative=1))
        zd = complex(mp.zeta(mp.mpc(0.5, t), derivative=1)); pred = -1j*np.exp(-1j*theta(np.array([t]))[0])*ex
        print(f"k={k+1} t={t:.6f} Z'_RS={zp:.6f} Z'_mp={ex:.6f} relerr={abs(zp-ex)/abs(ex):.2e} |zeta'-(-i e^-ith Z')|={abs(zd-pred):.1e}")
