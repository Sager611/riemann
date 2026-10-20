# pyright: reportMissingImports=false
# python-flint lives in a local venv (see README.md), not the system interpreter.

"""Galerkin discretization of the Weil functional on f supported in [-L/2, L/2].
Basis: orthonormal Legendre on [-a,a], a=L/2. Q(f) = poles + arch - primes (explicit formula,
h(r)=int g e^{irx}, g = f * f~). Verified conventions in check_ef.py."""
import sys, time, flint
from flint import arb, arb_mat, acb_mat, ctx

_gl = {}
def gauss(n):
    key = (n, ctx.prec)
    if key not in _gl:
        _gl[key] = [arb.legendre_p_root(n, k, weight=True) for k in range(n)]
    return _gl[key]

def leg_on(N, x):
    """sqrt((2j+1)/2) P_j(x), j<N (orthonormal on [-1,1])."""
    p0, p1 = arb(1), x
    out = [p0, p1][:N]
    for j in range(1, N-1):
        p0, p1 = p1, ((2*j+1)*x*p1 - j*p0)/(j+1)
        out.append(p1)
    return [out[j]*((arb(2*j+1)/2).sqrt()) for j in range(N)]

def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return arb(p).log() if m == 1 else None
    return None

def parseL(s):
    s = str(s)
    if s.startswith('log'): return arb(int(s[3:])).log()
    return arb(s)

def corr(N, a, t, Ny):
    """C_jk(t) = int phi_j(y) phi_k(y-t) dy, t in [0,2a]."""
    half = a - t/2
    if float(half) <= 0: return arb_mat(N, N)
    A = arb_mat(N, Ny); B = arb_mat(N, Ny)
    sc = 1/a
    for q, (x, w) in enumerate(gauss(Ny)):
        y = t/2 + half*x
        u = leg_on(N, y/a); v = leg_on(N, (y - t)/a)
        ww = w*half*sc
        for j in range(N):
            A[j, q] = u[j]*ww; B[j, q] = v[j]
    return A*B.transpose()

def weil_matrix(L, N, Q=None, verbose=False):
    a = L/2
    Q = Q or N + 120
    Ny = N + 1
    I = arb_mat(N, N)
    for j in range(N): I[j, j] = 1
    # archimedean
    M = I*(-(arb.const_euler() + arb.pi().log() + (1 - (-4*a).exp()).log()))
    for (x, w) in gauss(Q):
        t = a*(1 + x)
        C = corr(N, a, t, Ny)
        S = C + C.transpose()
        e2 = (-2*t).exp()
        M += (I*(2*e2) - S*((-t/2).exp())) * (w*a/(1 - e2))
    # poles
    Av = arb_mat(N, 1); Bv = arb_mat(N, 1)
    for (x, w) in gauss(N + 60):
        y = a*x
        u = leg_on(N, x)
        e = (y/2).exp()
        for j in range(N):
            Av[j, 0] += u[j]*e*w*a.sqrt()
    for j in range(N): Bv[j, 0] = Av[j, 0]*(-1)**j
    M += Av*Bv.transpose() + Bv*Av.transpose()
    # primes n < e^L
    Lf = float(L); n = 2
    import math
    while math.log(n) < Lf - 1e-15:
        lam = vm(n)
        if lam is not None:
            C = corr(N, a, arb(n).log(), Ny)
            M -= (C + C.transpose())*(lam/arb(n).sqrt())
        n += 1
    return M

def sub(M, idx):
    m = len(idx); R = arb_mat(m, m)
    for i, p in enumerate(idx):
        for j, q in enumerate(idx): R[i, j] = M[p, q].mid()
    return R

def lam_min(M):
    ev = acb_mat(M).eig(algorithm="approx")
    return min(e.real.mid() for e in ev)

if __name__ == "__main__":
    prec, N = int(sys.argv[1]), int(sys.argv[2])
    ctx.prec = prec
    for s in sys.argv[3:]:
        t0 = time.time()
        L = parseL(s)
        M = weil_matrix(L, N)
        le = lam_min(sub(M, list(range(0, N, 2)))); lo = lam_min(sub(M, list(range(1, N, 2))))
        print(f"L={s} ({float(L):.6f}) N={N} prec={prec} lam_even={le.str(8, radius=False)} lam_odd={lo.str(8, radius=False)} t={time.time()-t0:.1f}s", flush=True)
