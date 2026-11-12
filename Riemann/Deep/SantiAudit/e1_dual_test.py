# E1 numeric tests of the oscillatory dual term T(n) of (14).
# Test 1 (identity): Mellin form (14) vs x-space form -pi * int A(xi) xi^{-z/2} B(4 pi sqrt(n xi/(pq))) dxi,
#   B = sin(pi d/2) J_d + cos(pi d/2) Y_d, A = A_{-1,0}(xi) = xi^{-beta} exp(c^2/4), c = (1-2beta)/H + 2i(T-2 pi xi)/H.
# Test 2 (size): S = sum_n d(n)|T(n)| vs the E1 bound scale H (pq/T)^{1/4} N0^{3/4}, N0 = pq T/H^2,
#   and the per-(g,p,q) ratio r = (4 pi/(pq)) S / sqrt(T/H) (E1: r = O(T^eps)); decay of |T(n)| beyond N0.
import sys
mode = sys.argv[1]
if mode == "identity_mp":
    from mpmath import mp, mpf, mpc, quad, gamma, cos, sin, pi, exp, sqrt, besselj, bessely, linspace
    mp.dps = 20
    T, H = mpf(60), mpf(12)
    al, be = mpc('0.01', '0.02'), mpc('0.03', '-0.01'); z, d = al - be, al + be
    pq = 3
    A = lambda xi: xi**(-be) * exp(((1 - 2*be)/H + 2j*(T - 2*pi*xi)/H)**2 / 4)
    lo, hi = (T - 7*H)/(2*pi), (T + 7*H)/(2*pi)
    lo = max(lo, mpf('0.05'))
    nodes = linspace(lo, hi, 12)
    At = lambda s: quad(lambda xi: A(xi) * xi**(s - 1), nodes)
    sig = mpf('0.15')
    def TM(n):
        f = lambda t: (gamma(sig + 1j*t - d/2) * gamma(sig + 1j*t + d/2) * cos(pi*(sig + 1j*t))
                       * (mpf(pq)/(4*pi**2*n))**(sig + 1j*t) * At(1 - (sig + 1j*t) - z/2))
        return quad(f, linspace(-45, 45, 31)) / (2*pi)
    Bf = lambda x: sin(pi*d/2)*besselj(d, x) + cos(pi*d/2)*bessely(d, x)
    TX = lambda n: -pi * quad(lambda xi: A(xi) * xi**(-z/2) * Bf(4*pi*sqrt(n*xi/pq)), nodes)
    for n in (1, 2, 5):
        a, b = TM(n), TX(n)
        print("n", n, "Mellin", mp.nstr(a, 12), " x-space", mp.nstr(b, 12), " rel diff", mp.nstr(abs(a-b)/abs(b), 3))
else:
    import numpy as np
    from scipy.special import jv, yv
    al, be = 0.01, 0.02; z, d = al - be, al + be
    def divisors_count(N):
        dd = np.zeros(N + 1, dtype=np.int64)
        for k in range(1, N + 1): dd[k::k] += 1
        return dd
    print("T H pq N0 | S | S/pred | r=(4pi/pq)S/sqrt(T/H) | max_{n<=N0}|T|/(H X0^-1/2) | |T(n)| at n=4N0,8N0 rel. to n~N0")
    for (T, H, pq) in [(400., 40., 4), (400., 40., 16), (1600., 80., 16), (1600., 80., 64), (6400., 160., 64), (6400., 160., 256), (25600., 320., 256)]:
        N0 = pq*T/H**2; nmax = int(10*N0) + 20
        xi = np.linspace(max((T - 8*H)/(2*np.pi), 1e-3), (T + 8*H)/(2*np.pi), 3000); dx = xi[1] - xi[0]
        c = (1 - 2*be)/H + 2j*(T - 2*np.pi*xi)/H
        A = xi**(-be) * np.exp(c**2/4) * xi**(-z/2)
        dn = divisors_count(nmax)
        Tn = np.zeros(nmax + 1)
        for n in range(1, nmax + 1):
            X = 4*np.pi*np.sqrt(n*xi/pq)
            B = np.sin(np.pi*d/2)*jv(d, X) + np.cos(np.pi*d/2)*yv(d, X)
            Tn[n] = abs(-np.pi*np.sum(A*B)*dx)
        S = float(np.sum(dn[1:]*Tn[1:]))
        pred = H*(pq/T)**0.25*max(N0, 1.0)**0.75
        X0 = lambda n: 4*np.pi*np.sqrt(n*T/(2*np.pi*pq))
        m1 = max(Tn[n]/(H*X0(n)**-0.5) for n in range(1, max(2, int(N0)) + 1))
        ref = max(Tn[1:max(2, int(N0)) + 1])
        t4 = Tn[min(nmax, int(4*N0) + 1)]/ref; t8 = Tn[min(nmax, int(8*N0) + 1)]/ref
        print(int(T), int(H), pq, round(N0, 2), "| %.4g | %.3f | %.3f | %.3f | %.2e, %.2e" % (S, S/pred, 4*np.pi/pq*S/np.sqrt(T/H), m1, t4, t8))

if mode == "identity":
    # numpy/scipy version (real shifts): Mellin form (14) on Re u = 0.15 vs x-space Bessel form
    import numpy as np
    from scipy.special import gamma as G, jv, yv
    al, be = 0.01, 0.02; z, d = al - be, al + be
    for (T, H, pq) in [(60., 12., 3), (200., 25., 6)]:
        xi = np.linspace(max((T - 8*H)/(2*np.pi), 1e-3), (T + 8*H)/(2*np.pi), 4000); dx = xi[1] - xi[0]
        A = xi**(-be) * np.exp(((1 - 2*be)/H + 2j*(T - 2*np.pi*xi)/H)**2/4)
        sig = 0.15; Tmax = 14*T/H; t = np.linspace(-Tmax, Tmax, 40001); dt = t[1] - t[0]; u = sig + 1j*t
        At = np.empty(t.size, dtype=complex)
        for i0 in range(0, t.size, 400):
            uu = u[i0:i0+400][:, None]
            At[i0:i0+400] = np.sum(A[None, :] * xi[None, :]**(-uu - z/2), axis=1) * dx
        M = G(u - d/2) * G(u + d/2) * np.cos(np.pi*u)
        for n in (1, 2, 5, 9):
            TM = np.sum(M * (pq/(4*np.pi**2*n))**u * At) * dt / (2*np.pi)
            X = 4*np.pi*np.sqrt(n*xi/pq)
            TX = -np.pi*np.sum(A * xi**(-z/2) * (np.sin(np.pi*d/2)*jv(d, X) + np.cos(np.pi*d/2)*yv(d, X))) * dx
            print("T=%g H=%g pq=%d n=%d  Mellin %.10e%+.10ej  x-space %.10e%+.10ej  rel diff %.2e" % (T, H, pq, n, TM.real, TM.imag, TX.real, TX.imag, abs(TM-TX)/abs(TX)))
