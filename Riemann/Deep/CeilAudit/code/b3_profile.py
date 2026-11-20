# ceilaudit B3: SL/E1 dual term T(n) with the EXACT Gaussian profile (no Taylor in T/H^2), incl. theta < 1/2.
# T(n) = -pi int A(xi) xi^{-z/2} B(X(xi)) dxi, X = 4 pi sqrt(n xi/(pq)), A = xi^{-beta} pi^{-1/2} int e^{-y^2+(1-2beta)y/H} e^{2iwy} e^{-2 pi i xi g(y)} dy,
# g = e^{2y/H}-1-2y/H, w = (T-2 pi xi)/H  (= sum_j of E1Audit's A_j, eps = -1).  Substitution u = H(e^{2y/H}-1) makes the y-integral an FFT in u.
import sys, math, numpy as np
from scipy.special import roots_genlaguerre, gamma as cgamma
NL = 80
def apm(X, d):
    t, wt = roots_genlaguerre(NL, d.real - 0.5)
    base = wt * np.exp(1j * d.imag * np.log(t)) / cgamma(d + 0.5)
    tt = t[None, :] / (2 * X[:, None])
    return ((1 + 1j * tt) ** (d - 0.5)) @ base, ((1 - 1j * tt) ** (d - 0.5)) @ base
def Bfun(X, d):
    out = np.empty(len(X), complex)
    for i in range(0, len(X), 4000):
        x = X[i:i+4000]; ap, am = apm(x, d)
        out[i:i+4000] = np.sqrt(2 / (np.pi * x)) * (-0.5j) * (np.exp(1j * (x - np.pi / 4)) * ap - np.exp(-1j * (x - np.pi / 4)) * am)
    return out
def profile_exact(T, H, beta, du, N):
    # returns w-grid and Phi(w) = pi^{-1/2} int e^{-y^2+(1-2b)y/H} e^{2iwy} e^{-i(T-Hw) g(y)} dy
    m = np.arange(N); u0 = -N * du / 2; u = u0 + m * du
    ok = u > -H * 0.999
    y = np.where(ok, 0.5 * H * np.log1p(np.where(ok, u, 0) / H), -1e9)
    amp = np.where(np.abs(y) < 9, np.exp(-y ** 2 + (1 - 2 * beta) * y / H), 0) * np.where(ok, 0.5 / (1 + np.where(ok, u, 0) / H), 0)
    gph = np.where(ok, u / H - np.log1p(np.where(ok, u, 0) / H), 0)   # g(y) = u/H - log(1+u/H)
    G = amp * np.exp(-1j * T * gph)
    k = np.fft.fftfreq(N) * N; w = 2 * np.pi * k / (N * du)
    Phi = du * np.exp(1j * w * u0) * N * np.fft.ifft(G) / math.sqrt(math.pi)
    o = np.argsort(w); return w[o], Phi[o]
def profile_taylor(T, H, beta, w, J):
    y = np.arange(-12, 12.0001, 0.02); g = np.expm1(2 * y / H) - 2 * y / H; xi = (T - H * w) / (2 * np.pi)
    out = np.zeros(len(w), complex)
    for j in range(J):
        amp = np.exp(-y ** 2 + (1 - 2 * beta) * y / H) * g ** j * 0.02 / (math.factorial(j) * math.sqrt(math.pi))
        out += (-2j * np.pi * xi) ** j * (np.exp(2j * np.outer(w, y)) @ amp)
    return out
def brute(T, H, beta, w, dy=2e-4):
    y = np.arange(-9, 9, dy); g = np.expm1(2 * y / H) - 2 * y / H
    return np.array([np.sum(np.exp(-y ** 2 + (1 - 2 * beta) * y / H + 2j * wi * y - 1j * (T - H * wi) * g)) * dy / math.sqrt(math.pi) for wi in w])
def setup(T, H, beta, alpha):
    tau = T / H ** 2
    umax = H * (math.exp(18 / H) - 1)
    du = min(2e-3, 2 * np.pi / (8 * max(tau, 0.05) * umax))
    N = 1 << int(math.ceil(math.log2(2 * np.pi / (0.02 * du))))   # dw = 2pi/(N du) <= 0.02
    w, Phi = profile_exact(T, H, beta, du, N)
    Wt = min(8 * math.sqrt(1 + 4 * tau ** 2) + 12, 0.5 * T / H)
    keep = (np.abs(w) <= Wt)
    return w[keep].copy(), Phi[keep].copy(), N, du
def Tn(T, H, pq, ns, alpha, beta, w, Phi):
    z, d = alpha - beta, alpha + beta
    dw = w[1] - w[0]; xi = (T - H * w) / (2 * np.pi)
    A = xi ** (-beta) * Phi * xi ** (-z / 2)
    return np.array([-np.pi * (H / (2 * np.pi)) * dw * np.sum(A * Bfun(4 * np.pi * np.sqrt(n * xi / pq), d)) for n in ns])
def divcount(n):
    return sum(2 if k * k != n else 1 for k in range(1, int(math.isqrt(n)) + 1) if n % k == 0)
if __name__ == '__main__':
    al, be = complex(sys.argv[1]) if len(sys.argv) > 1 else 0.01, complex(sys.argv[2]) if len(sys.argv) > 2 else 0.02
    al, be = complex(al), complex(be)
    print(f"shifts alpha={al}, beta={be}", flush=True)
    # validation 1: exact FFT profile vs Taylor (theta > 1/2) and vs brute y-quadrature (theta < 1/2)
    for T, H in [(1e6, 3000.0), (1e6, 250.0)]:
        w, Phi, N, du = setup(T, H, be.real, al)
        ws = np.array([0.0, 1.3, -2.7, 5.1, 0.5 * w[-1]])
        ex = np.interp(ws, w, Phi.real) + 1j * np.interp(ws, w, Phi.imag)
        bf = brute(T, H, be.real, ws)
        msg = f"validate T={T:.0e} H={H} tau={T/H**2:.3g} N={N} du={du:.2e} |w|<={w[-1]:.1f}: max|FFT-brute|/max|Phi| = {np.abs(ex-bf).max()/np.abs(Phi).max():.1e}"
        if T / H ** 2 < 0.5:
            ta = profile_taylor(T, H, be.real, ws, 8); msg += f"; |FFT-Taylor(J=8)|/max = {np.abs(ex-ta).max()/np.abs(Phi).max():.1e}"
        print(msg, flush=True)
    rows = [tuple(map(float, r.split(':'))) for r in sys.argv[3].split(',')] if len(sys.argv) > 3 else [(1e6, 1000.0), (1e6, 500.0), (1e6, 250.0), (1e6, 125.0), (1e8, 10000.0), (1e8, 2500.0), (1e8, 1250.0)]
    for T, H in rows:
        w, Phi, N, du = setup(T, H, be.real, al)
        tau = T / H ** 2; L1 = (H / (2 * np.pi)) * (w[1] - w[0]) * np.abs(Phi).sum()
        for pq in [6, 30]:
            N0 = pq * T / H ** 2
            full = 12 * N0 <= 300
            ns = list(range(1, int(math.ceil(12 * N0)) + 1)) if full else sorted(set([1, 2, 3] + [max(1, int(round(N0 * f))) for f in (0.25, 0.5, 1, 2, 4, 8)]))
            v = Tn(T, H, pq, ns, al, be, w, Phi)
            bd = np.array([H * (2 * math.sqrt(2 * math.pi) * math.sqrt(n * T / pq)) ** -0.5 for n in ns])
            r = np.abs(v) / bd
            mx = max(r[i] for i, n in enumerate(ns) if n <= max(1, N0))
            dec = {f: max([r[i] for i, n in enumerate(ns) if n >= f * N0] or [0]) / mx for f in (2, 4, 8)}
            S = (sum(divcount(n) * abs(v[i]) for i, n in enumerate(ns)) / (H * (pq / T) ** 0.25 * max(N0, 1) ** 0.75)) if full else float('nan')
            print(f"T={T:.0e} H={H:.0f} theta={math.log(H)/math.log(T):.3f} tau=T/H^2={tau:.3g} pq={pq} N0={N0:.3g}: L1(A)/H={L1/H:.3f} (sqrt(tau)={math.sqrt(tau):.2f}); "
                  f"max_(n<=N0)|T(n)|/(H X0^-1/2)={mx:.3f}; decay >=2N0 {dec[2]:.1e} >=4N0 {dec[4]:.1e} >=8N0 {dec[8]:.1e}; S/pred={S:.3f}", flush=True)
