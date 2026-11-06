# zeros.py (explore23): zeta zeros 0 < gamma <= TMAX by numpy Riemann-Siegel (C0..C4), validated against mpmath.
import sys, time, json, numpy as np, mpmath as mp
mp.mp.dps = 50
PI = np.pi
def Psi(p):
    return mp.cos(2*mp.pi*(p*p - p - mp.mpf(1)/16)) / mp.cos(2*mp.pi*p)
# Taylor coefficients of Psi around p = 1/2 (entire); Psi has a removable 0/0 at p=1/4,3/4 only.
NT = 44
c = mp.taylor(Psi, mp.mpf('0.5'), NT)
c = np.array([float(x) for x in c])
fact = np.array([float(mp.factorial(j)) for j in range(NT + 1)])
def dPsi(k, p):
    # k-th derivative of Psi at p (array), via Taylor series around 1/2
    x = p - 0.5
    s = np.zeros_like(p)
    for j in range(NT, k - 1, -1):
        s = s * x + c[j] * fact[j] / fact[j - k]
    return s
def theta(t):
    t = np.asarray(t, dtype=float)
    return t/2*np.log(t/(2*PI)) - t/2 - PI/8 + 1/(48*t) + 7/(5760*t**3) + 31/(80640*t**5)
def Z(t):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.empty_like(t)
    a = np.sqrt(t/(2*PI)); N = np.floor(a).astype(int); p = a - N
    th = theta(t)
    Nmax = N.max()
    n = np.arange(1, Nmax + 1, dtype=float)
    CH = 4096
    for i in range(0, len(t), CH):
        tt = t[i:i+CH]; NN = N[i:i+CH]
        M = n[None, :] <= NN[:, None]
        terms = np.cos(th[i:i+CH, None] - tt[:, None]*np.log(n)[None, :]) / np.sqrt(n)[None, :]
        out[i:i+CH] = 2*np.sum(np.where(M, terms, 0.0), axis=1)
    pi2 = PI*PI
    C0 = dPsi(0, p)
    C1 = -dPsi(3, p)/(96*pi2)
    C2 = dPsi(2, p)/(64*pi2) + dPsi(6, p)/(18432*pi2*pi2)
    C3 = -dPsi(1, p)/(64*pi2) - dPsi(5, p)/(3840*pi2*pi2) - dPsi(9, p)/(5308416*pi2**3)
    C4 = dPsi(0, p)/(128*pi2) + 19*dPsi(4, p)/(24576*pi2*pi2) + 11*dPsi(8, p)/(5898240*pi2**3) + dPsi(12, p)/(2038431744*pi2**4)
    u = 1/a
    R = (-1.0)**(N - 1) * np.sqrt(u) * (C0 + C1*u + C2*u**2 + C3*u**3 + C4*u**4)
    return out + R
if __name__ == '__main__':
    TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 20000.0
    OUT = sys.argv[2] if len(sys.argv) > 2 else 'zeta_zeros.npy'
    t0 = time.time()
    mp.mp.dps = 20
    for tv in [60.0, 200.0, 1000.0, 7005.06, 10000.0, 19999.3]:
        print('check Z', tv, float(Z(tv)[0]), float(mp.siegelz(tv)), flush=True)
    # grid
    T0 = 10.0
    ts = [T0]
    t = T0
    while t < TMAX:
        sp = 2*PI/np.log(max(t, 20)/(2*PI))
        t += sp/20
        ts.append(t)
    ts = np.array(ts)
    zv = Z(ts)
    print('grid', len(ts), round(time.time()-t0, 1), 's', flush=True)
    idx = np.nonzero(np.sign(zv[:-1]) * np.sign(zv[1:]) < 0)[0]
    lo = ts[idx].copy(); hi = ts[idx+1].copy(); flo = zv[idx].copy()
    for it in range(60):
        mid = (lo + hi)/2; fm = Z(mid)
        s = np.sign(fm) == np.sign(flo)
        lo = np.where(s, mid, lo); flo = np.where(s, fm, flo); hi = np.where(s, hi, mid)
    zs = (lo + hi)/2
    print('found', len(zs), round(time.time()-t0, 1), 's', flush=True)
    # Turing-type count check: theta(g_n)/pi + 1 - (n - 1/2) = -S(g_n) approx; a missed pair shifts by 2
    nn = np.arange(1, len(zs) + 1)
    dev = theta(zs)/PI + 1 - (nn - 0.5)
    print('dev min/max', dev.min(), dev.max(), flush=True)
    # local search for missed pairs: windows where |Z| has a small local extremum without sign change
    absz = np.abs(zv)
    lm = np.nonzero((absz[1:-1] < absz[:-2]) & (absz[1:-1] < absz[2:]))[0] + 1
    sus = []
    for j in lm:
        if np.sign(zv[j-1]) == np.sign(zv[j+1]) == np.sign(zv[j]):
            # fine scan
            f = np.linspace(ts[j-1], ts[j+1], 401); fz = Z(f)
            ch = np.nonzero(np.sign(fz[:-1])*np.sign(fz[1:]) < 0)[0]
            if len(ch):
                sus.append((ts[j], len(ch)))
                for k in ch:
                    a_, b_ = f[k], f[k+1]; fa = fz[k]
                    for it in range(60):
                        m_ = (a_+b_)/2; fmm = Z(m_)[0]
                        if np.sign(fmm) == np.sign(fa): a_, fa = m_, fmm
                        else: b_ = m_
                    zs = np.append(zs, (a_+b_)/2)
    zs = np.sort(zs)
    print('missed found by fine scan', sus[:20], flush=True)
    nn = np.arange(1, len(zs) + 1)
    dev = theta(zs)/PI + 1 - (nn - 0.5)
    print('after: n', len(zs), 'dev min/max', round(dev.min(),4), round(dev.max(),4), 'mean', round(dev.mean(),5), flush=True)
    np.save(OUT, zs)
    # validate vs mpmath zetazero
    mp.mp.dps = 20
    for k in [1, 2, 100, 1000, 6709, 6710, 10000]:
        if k <= len(zs):
            g = mp.zetazero(k).imag
            print('zetazero', k, float(g), zs[k-1], float(g) - zs[k-1], flush=True)
    print('done', round(time.time()-t0, 1), 's', flush=True)
