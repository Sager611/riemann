# dh.py (explore23): Davenport-Heilbronn zeros up to TMAX: on-line by sign changes of Z_DH, off-line located from the count deficit.
# DH: Lambda(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s), f = 5^{-s} sum_r a_r zeta(s, r/5), a = (1, k, -k, -1), k = kappa (as explore22 ecommon.py).
import sys, json, time, numpy as np, mpmath as mp
PI = np.pi
k = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1)
A = {1: 1.0, 2: k, 3: -k, 4: -1.0}
B2 = [mp.bernoulli(2*j) for j in range(1, 13)]
B2 = np.array([float(b/mp.factorial(2*j)) for j, b in enumerate(B2, start=1)])
def hurwitz(s, a, M):
    # Euler-Maclaurin, vector s (complex), scalar a
    kk = np.arange(M, dtype=float) + a
    out = np.zeros(len(s), dtype=complex)
    for i in range(0, len(s), 2048):
        ss = s[i:i+2048]
        out[i:i+2048] = np.exp(-ss[:, None]*np.log(kk)[None, :]).sum(1)
    N = M + a
    out += N**(1 - s)/(s - 1) + 0.5*N**(-s)
    fac = s.copy(); term = N**(-s - 1)
    for j in range(len(B2)):
        out += B2[j]*fac*term
        fac = fac*(s + 2*j + 1)*(s + 2*j + 2)
        term = term/(N*N)
    return out
def f(s, M=None):
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    if M is None: M = int(np.max(np.abs(s.imag))/2) + 40
    return 5.0**(-s)*sum(A[r]*hurwitz(s, r/5.0, M) for r in A)
def vartheta(t):
    return np.array([float(mp.im(mp.loggamma(mp.mpf(0.75) + 0.5j*mp.mpf(x)))) + 0.5*x*np.log(5/PI) for x in t])
def Zdh(t):
    t = np.atleast_1d(t)
    return (np.exp(1j*vartheta(t))*f(0.5 + 1j*t)).real
if __name__ == '__main__':
    TMAX = float(sys.argv[1]); t0 = time.time()
    mp.mp.dps = 20
    kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
    def fmp(s):
        return mp.power(5, -s)*(mp.zeta(s, mp.mpf(1)/5) + kap*mp.zeta(s, mp.mpf(2)/5) - kap*mp.zeta(s, mp.mpf(3)/5) - mp.zeta(s, mp.mpf(4)/5))
    for s in [0.5 + 20j, 0.5 + 150.3j, 0.8 + 300.1j, 0.5 + TMAX*1j - 1j]:
        a = f(np.array([s]))[0]; b = complex(fmp(mp.mpc(s.real, s.imag)))
        print('check f', s, abs(a - b), abs(b), flush=True)
    ts = np.arange(1.0, TMAX, 0.02)
    z = Zdh(ts)
    print('grid', len(ts), round(time.time() - t0, 1), flush=True)
    idx = np.nonzero(np.sign(z[:-1])*np.sign(z[1:]) < 0)[0]
    lo = ts[idx]; hi = ts[idx + 1]; flo = z[idx]
    for it in range(45):
        mid = (lo + hi)/2; fm = Zdh(mid)
        sm = np.sign(fm) == np.sign(flo)
        lo = np.where(sm, mid, lo); flo = np.where(sm, fm, flo); hi = np.where(sm, hi, mid)
    on = (lo + hi)/2
    vt = vartheta(on)
    n = np.arange(1, len(on) + 1)
    dev = vt/PI - (n - 0.5)
    print('on-line', len(on), 'dev first/last', dev[:3], dev[-3:], round(time.time() - t0, 1), flush=True)
    # off-line search: vectorized numpy Newton from a dense grid of starts (sigma 0.55..1.9, t step 0.25), then mpmath polish
    T_, S_ = np.meshgrid(np.arange(2.0, TMAX, 0.25), np.array([0.56, 0.65, 0.8, 1.0, 1.25, 1.6, 1.9]))
    s = (S_ + 1j*T_).ravel()
    for it in range(40):
        h = 1e-6
        fs = f(s, M=int(TMAX/2) + 60); d = (f(s + h, M=int(TMAX/2) + 60) - f(s - h, M=int(TMAX/2) + 60))/(2*h)
        step = fs/d
        step = np.where(np.abs(step) > 0.5, 0.5*step/np.abs(step), step)
        s = s - step
        s = np.where((s.real < 0.5) | (s.real > 3.0) | (s.imag < 0) | (s.imag > TMAX + 5), np.nan + 0j, s)
    ok = np.isfinite(s) & (np.abs(f(np.nan_to_num(s), M=int(TMAX/2) + 60)) < 1e-9) & (s.real > 0.5 + 1e-5)
    cand = s[ok]
    roots = []
    for r0 in cand:
        if all(abs(r0 - complex(q)) > 1e-5 for q in roots):
            r = mp.findroot(fmp, mp.mpc(r0.real, r0.imag), tol=1e-20)
            if all(abs(r - q) > 1e-8 for q in roots):
                roots.append(r)
    roots = sorted(roots, key=lambda r: float(r.imag))
    print('off-line (Re>1/2):', [(round(float(r.real), 6), round(float(r.imag), 6)) for r in roots], round(time.time() - t0, 1), flush=True)
    json.dump({'TMAX': TMAX, 'online': [float(x) for x in on], 'offline': [[float(r.real), float(r.imag)] for r in roots],
               'dev_head': [float(x) for x in dev[:5]], 'dev_tail': [float(x) for x in dev[-5:]]}, open('dh_zeros.json', 'w'))
    print('done', round(time.time() - t0, 1), flush=True)
