# Part 2b: I = (1/T) int_T^{2T} |zeta M(sigma+it)|^2 dt for Baluyot's mollifier, vs M1 + M2 and vs 1.
# usage: meanvalue.py T thetas(comma) h eps [extra sigmas comma]
import sys, time, json
import numpy as np
from zeta_rs import zeta_rs
from mollifier import coeffs, main_terms, m2_from_s2

T = float(sys.argv[1])
thetas = [float(v) for v in sys.argv[2].split(',')]
h = float(sys.argv[3])
eps = float(sys.argv[4])
extra = [float(v) for v in sys.argv[5].split(',')] if len(sys.argv) > 5 else [0.55, 0.6, 0.7, 0.8]
sig = set(extra)
for th in thetas:
    L = th * np.log(T)
    for c in (1, 2):
        s = 0.5 + c / L
        if s < 0.95:
            sig.add(round(s, 6))
sigmas = sorted(sig)
J = int(round(T / h))
t = T + h * np.arange(J + 1)
w = np.full(J + 1, h)
w[0] = w[-1] = h / 2
res = []
for sigma in sigmas:
    t0 = time.time()
    z = np.empty(J + 1, dtype=complex)
    for i in range(0, J + 1, 256):
        z[i:i + 256] = zeta_rs(sigma, t[i:i + 256])
    tz = time.time() - t0
    for th in thetas:
        t1 = time.time()
        y = T ** th
        N = max(1, int(np.floor(y ** (1 + eps))))
        a, mu, primes, W = coeffs(y, sigma, N)
        n = np.arange(1, N + 1, dtype=float)
        logn = np.log(n)
        an = a[1:] * n ** (-sigma)
        sel = an != 0
        an, logn = an[sel], logn[sel]
        B = int(max(8, min(128, 2_000_000 // max(1, len(logn)))))  # E is B x len(logn) complex: <= 32 MB
        kk = np.arange(B)
        E = np.exp(-1j * h * np.outer(kk, logn))
        I = 0.0
        for i in range(0, J + 1, B):
            m = min(B, J + 1 - i)
            v = an * np.exp(-1j * t[i] * logn)
            Mb = E[:m] @ v
            I += float((w[i:i + m] * np.abs(z[i:i + m] * Mb) ** 2).sum())
        I /= T
        M1, S2 = main_terms(a, mu, primes, sigma)
        M2 = m2_from_s2(S2, sigma, T)
        yp = y ** (1 - 2 * sigma)
        inrange = sigma >= 0.5 + 1 / np.log(y)
        r = dict(T=T, theta=th, y=y, N=N, sigma=sigma, c=(sigma - 0.5) * np.log(y), inrange=bool(inrange),
                 I=I, M1=M1, M2=M2, yp=yp, r1=(I - 1) / yp, off=I - M1 - M2,
                 bound=T ** -0.5 * y ** 0.875 + T ** -1 * y ** 1.75, h=h, eps=eps,
                 sec_zeta=tz, sec=time.time() - t1)
        res.append(r)
        print("T=%g th=%.2f y=%.1f N=%d sig=%.4f c=%.2f in=%d  I=%.6f M1=%.6f M2=%+.6f  (I-1)/yp=%.3f  I-M1-M2=%+.2e  [T^-1/2 y^7/8+T^-1 y^7/4=%.2e]  %.0fs"
              % (T, th, y, N, sigma, r['c'], inrange, I, M1, M2, r['r1'], r['off'], r['bound'], r['sec'] + tz), flush=True)
json.dump(res, open('../logs/meanvalue_T%g_h%g.json' % (T, h), 'w'), indent=1)
