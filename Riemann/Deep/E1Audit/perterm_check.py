# E1Audit: per-term test of E1 Steps 3-4: |T_j(n)| << (T/H^2)^j H X0^{-1/2} (1+sqrt(n/N0))^{-K}, complex shifts, j=0,1,2, T up to 1e7.
# T_j(n) = -pi int A_j(xi) xi^{-z/2} B(X(xi)) dxi, A_j = xi^{-beta} (-2 pi i xi)^j Phi_j((T-2 pi xi)/H), X = 4 pi sqrt(n xi/(pq)).
# Phi_j by trapezoid in y (expm1 bracket, spectrally accurate); B by the exact Hankel form with a_pm by generalized Gauss-Laguerre
# (validated against mpmath Bessel J_d, Y_d first). Trapezoid in w = (T-2 pi xi)/H on [-10,10].
import sys, math, numpy as np
from scipy.special import roots_genlaguerre, gamma as cgamma
NL = 80
def apm(X, d):
    t, wt = roots_genlaguerre(NL, d.real - 0.5)
    base = wt * np.exp(1j * d.imag * np.log(t)) / cgamma(d + 0.5)
    tt = t[None, :] / (2 * X[:, None])
    ap = ((1 + 1j * tt) ** (d - 0.5)) @ base; am = ((1 - 1j * tt) ** (d - 0.5)) @ base
    return ap, am
def Bfun(X, d):
    ap, am = apm(X, d)
    return np.sqrt(2 / (np.pi * X)) * (-0.5j) * (np.exp(1j * (X - np.pi / 4)) * ap - np.exp(-1j * (X - np.pi / 4)) * am)
def validate():
    import mpmath as m
    out = []
    for d in [complex(0.02, 0.17), 0j, complex(-0.2, 0.0), complex(0.1, -0.2)]:
        for X in [20.0, 100.0, 1000.0]:
            B = Bfun(np.array([X]), d)[0]; dm = m.mpc(d.real, d.imag)
            Bm = complex(m.sin(m.pi * dm / 2) * m.besselj(dm, X) + m.cos(m.pi * dm / 2) * m.bessely(dm, X))
            out.append(abs(B - Bm) / abs(Bm))
    print(f"validate Hankel-Laguerre B vs mpmath (d in {{0.02+0.17i,0,-0.2,0.1-0.2i}}, X in {{20,100,1000}}): max rel err {max(out):.2e}", flush=True)
def Phi(w, j, beta, H):
    y = np.arange(-12, 12.0001, 0.04); g = np.expm1(2 * y / H) - 2 * y / H
    amp = np.exp(-y ** 2 + (1 - 2 * beta) * y / H) * g ** j * 0.04 / (math.factorial(j) * math.sqrt(math.pi))
    out = np.empty(len(w), complex)
    for i in range(0, len(w), 500):
        out[i:i + 500] = np.exp(2j * np.outer(w[i:i + 500], y)) @ amp
    return out
def Tj(T, H, pq, nlist, j, alpha, beta):
    z, d = alpha - beta, alpha + beta
    w = np.linspace(-10, 10, 4001); dw = w[1] - w[0]; xi = (T - H * w) / (2 * np.pi)
    A = xi ** (-beta) * (-2j * np.pi * xi) ** j * Phi(w, j, beta, H) * xi ** (-z / 2)
    res = []
    for n in nlist:
        X = 4 * np.pi * np.sqrt(n * xi / pq)
        res.append(-np.pi * (H / (2 * np.pi)) * dw * np.sum(A * Bfun(X, d)))
    return np.array(res)
def divcount(n):
    return sum(2 if k * k != n else 1 for k in range(1, int(math.isqrt(n)) + 1) if n % k == 0)
def row(T, th, pq, j, alpha, beta, H=None, full=True):
    H = H or T ** th; N0 = pq * T / H ** 2
    nmax = max(4, int(math.ceil(12 * N0)))
    ns = list(range(1, nmax + 1)) if (full and nmax <= 2500) else sorted(set([1, 2, 3] + [max(1, int(round(N0 * f))) for f in (0.25, 0.5, 1, 2, 4, 8, 12)]))
    v = Tj(T, H, pq, ns, j, alpha, beta)
    bd = lambda n: (T / H ** 2) ** j * H * (2 * math.sqrt(2 * math.pi) * math.sqrt(n * T / pq)) ** -0.5
    r = np.array([abs(v[i]) / bd(n) for i, n in enumerate(ns)])
    below = [r[i] for i, n in enumerate(ns) if n <= max(1, N0)]
    mx = max(below)
    dec = {f: max([r[i] for i, n in enumerate(ns) if n >= f * N0] or [0]) / mx for f in (2, 4, 8)}
    S = sum(divcount(n) * abs(v[i]) for i, n in enumerate(ns)) if ns == list(range(1, nmax + 1)) else float('nan')
    pred = (T / H ** 2) ** j * H * (pq / T) ** 0.25 * max(N0, 1) ** 0.75 if N0 >= 1 else (T / H ** 2) ** j * H * (pq / T) ** 0.25
    print(f"T={T:.0e} th={th} H={H:.1f} pq={pq} N0={N0:.3g} j={j}: max_(n<=N0) |T|/bound={mx:.3f}; decay rel. max: >=2N0 {dec[2]:.1e}, >=4N0 {dec[4]:.1e}, >=8N0 {dec[8]:.1e}; "
          f"S=sum d(n)|T|/[(T/H^2)^j H (pq/T)^1/4 N0^3/4]={S / pred:.3f}", flush=True)
if __name__ == '__main__':
    validate()
    print("cross-check vs E1 e1_dual_size.log (real shifts a=0.01, b=0.02, j=0): expected 0.107 (T=400,pq=4) and 0.545 (T=25600,H=320,pq=256)")
    row(400, 0, 4, 0, 0.01, 0.02, H=40.0); row(25600, 0, 256, 0, 0.01, 0.02, H=320.0)
    al, be = complex(0.05, 0.1), complex(-0.03, 0.07)
    print(f"complex shifts alpha={al}, beta={be} (d={al+be}, z={al-be})")
    for T in [1e5, 1e6, 1e7]:
        for th in [0.505, 0.55]:
            for pq in [6, 210]:
                for j in [0, 1, 2]:
                    row(T, th, pq, j, al, be)
    print("d = 0 exactly (alpha=0.01+0.02i, beta=-0.01-0.02i), j=0,1:")
    for T, th, pq in [(1e6, 0.505, 30), (1e7, 0.55, 210)]:
        for j in [0, 1]: row(T, th, pq, j, complex(0.01, 0.02), complex(-0.01, -0.02))
    print("N0 < 1 regime (p=q=1):")
    for T, th in [(1e6, 0.55), (1e6, 0.6), (1e7, 0.75)]:
        H = T ** th; N0 = T / H ** 2; v = Tj(T, H, 1, [1, 2, 5], 0, al, be)
        bd = H * (2 * math.sqrt(2 * math.pi) * math.sqrt(T)) ** -0.5
        print(f"  T={T:.0e} th={th} N0={N0:.3g}: |T(1)|,|T(2)|,|T(5)| / (H X0(1)^-1/2) = {abs(v[0])/bd:.2e}, {abs(v[1])/bd:.2e}, {abs(v[2])/bd:.2e}")
