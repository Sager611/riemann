# dh_ft.py (explore35): Davenport-Heilbronn on-line zeros up to TMAX (sign changes of Z_DH, ExploreHelio's Euler-Maclaurin f),
# then sum_on x^{i gamma} vs DH's Landau term -(T/2pi) c(x)/sqrt(x), c = coefficients of -f'/f (c * a = a log), and the off-line deficit.
import sys, json, time, math, numpy as np, mpmath as mp
PI = np.pi
k = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1)
A = {1: 1.0, 2: k, 3: -k, 4: -1.0}
B2 = np.array([float(mp.bernoulli(2*j)/mp.factorial(2*j)) for j in range(1, 13)])
def hurwitz(s, a, M, CH=512):
    kk = np.log(np.arange(M, dtype=float) + a)
    out = np.zeros(len(s), dtype=complex)
    for i in range(0, len(s), CH):
        ss = s[i:i+CH]
        out[i:i+CH] = np.exp(-ss[:, None]*kk[None, :]).sum(1)
    N = M + a
    out += N**(1 - s)/(s - 1) + 0.5*N**(-s)
    fac = s.copy(); term = N**(-s - 1)
    for j in range(len(B2)):
        out += B2[j]*fac*term; fac = fac*(s + 2*j + 1)*(s + 2*j + 2); term = term/(N*N)
    return out
def f(s, M):
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    return 5.0**(-s)*sum(A[r]*hurwitz(s, r/5.0, M) for r in A)
mp.mp.dps = 20
def vartheta(t):
    return np.array([float(mp.im(mp.loggamma(mp.mpf(0.75) + 0.5j*mp.mpf(float(x))))) + 0.5*float(x)*np.log(5/PI) for x in t])
def Zdh(t, M):
    t = np.atleast_1d(t)
    return (np.exp(1j*vartheta(t))*f(0.5 + 1j*t, M)).real
def coeffs(nmax):
    a = [0.0]*(nmax + 1)
    for n in range(1, nmax + 1): a[n] = [0.0, 1.0, k, -k, -1.0][n % 5]
    c = [0.0]*(nmax + 1)
    for n in range(2, nmax + 1):
        c[n] = a[n]*math.log(n) - sum(c[d]*a[n//d] for d in range(2, n) if n % d == 0)
    return a, c
if __name__ == '__main__':
    TMAX = float(sys.argv[1]); t0 = time.time(); M = int(TMAX/2) + 60
    ts = np.arange(1.0, TMAX, 0.02); z = np.empty(len(ts))
    for i in range(0, len(ts), 20000): z[i:i+20000] = Zdh(ts[i:i+20000], M)
    print('grid', len(ts), round(time.time() - t0, 1), flush=True)
    idx = np.nonzero(np.sign(z[:-1])*np.sign(z[1:]) < 0)[0]
    lo = ts[idx]; hi = ts[idx + 1]; flo = z[idx]
    for it in range(40):
        mid = (lo + hi)/2; fm = Zdh(mid, M); sm = np.sign(fm) == np.sign(flo)
        lo = np.where(sm, mid, lo); flo = np.where(sm, fm, flo); hi = np.where(sm, hi, mid)
    on = (lo + hi)/2
    vt = vartheta(on); n = np.arange(1, len(on) + 1); dev = vt/PI - (n - 0.5)
    print('on-line', len(on), 'dev at T: %.3f' % dev[-1], round(time.time() - t0, 1), flush=True)
    # deficit profile: dev ~ S-type fluctuation + (number of zeros off the line below t)
    marks = {}
    for T in (420, 1000, 1500, 2000, 2500, 3000):
        if T <= TMAX:
            j = np.searchsorted(on, T) - 1; marks[T] = float(dev[j])
    print('dev profile', {kk: round(v, 2) for kk, v in marks.items()}, flush=True)
    a, c = coeffs(40)
    rows = []
    for x in range(2, 41):
        ph = np.mod(on*math.log(x), 2*PI); S = complex(np.cos(ph).sum(), np.sin(ph).sum())
        P = -(TMAX/(2*PI))*c[x]/math.sqrt(x)
        rows.append((x, a[x], c[x], S.real, S.imag, P, S.real - P))
    print('| x | a(x) | c(x) | Re S_on | Im S_on | Landau -(T/2pi) c/sqrt(x) | Re S_on - P |')
    for r in rows: print('| %d | %+.4f | %+.4f | %+.1f | %+.1f | %+.1f | %+.1f |' % r, flush=True)
    json.dump({'TMAX': TMAX, 'online': [float(x) for x in on], 'dev_profile': marks, 'rows': rows}, open(sys.argv[2], 'w'))
    print('done', round(time.time() - t0, 1), flush=True)
