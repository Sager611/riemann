# JWAudit spot check of one certified SavvaK cell by a different method (no saddle expansion, no Q_k, no R_k series):
# true h(n0+j) = log Mt(n0+j) by direct quadrature (real axis, k = 1 theta term; k >= 2 is < e^{-10^38} relative),
# exact six-match system Delta^2 h(n0+k) = f(A+k)+f(C+k)-f(B+k)+f(n0+1/2+k)-f(D+k) (k = 0..3), f(U) = log(1+1/U),
# solved by Newton at 330 digits; then n0^5 L E_F^(6) from Delta^6 h and mpmath psi^(5).
# Run: uv run --no-project --with mpmath python3 spotcheck.py [log10 n0]
import mpmath as mp, sys, time
mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 330
LS = int(sys.argv[3]) if len(sys.argv) > 3 else 80; NE = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
pi = mp.pi
e10 = int(sys.argv[1]) if len(sys.argv) > 1 else 40
n0 = mp.mpf(10) ** e10
t0 = time.time()
def Z(u): return u * (pi * mp.exp(2 * u) - mp.mpf(9) / 4)
def Zp(u): return pi * mp.exp(2 * u) * (1 + 2 * u) - mp.mpf(9) / 4
def u0of(x):
    u = mp.log(x) / 2
    for _ in range(600):
        du = (Z(u) - x) / Zp(u); u -= du
        if abs(du) < mp.mpf(10) ** (-mp.mp.dps + 12): break
    return u
def f(u, x): return 2 * x * mp.log(u) + mp.mpf(9) * u / 2 - pi * mp.exp(2 * u)
def f2(u, x): return -2 * x / u ** 2 - 4 * pi * mp.exp(2 * u)
def amp(u): return 1 - 3 / (2 * pi) * mp.exp(-2 * u)
def h(x):   # log Mt(x) - log(2 pi^2)
    u0 = u0of(x); r = mp.sqrt(-f2(u0, x)); F0 = f(u0, x)
    g = lambda sg: mp.exp(f(u0 + sg / r, x) - F0) * amp(u0 + sg / r)
    I, err = mp.quad(g, [-46, -30, -20, -10, -4, 0, 4, 10, 20, 30, 46], error=True)
    return F0 + mp.log(I / r), err / abs(I)
H = []; qerr = 0
for j in range(7):
    v, e = h(n0 + j); H.append(v); qerr = max(qerr, e)
print("quadrature: max relative error estimate %s, time %.1f s" % (mp.nstr(qerr, 3), time.time() - t0), flush=True)
L = 2 * u0of(n0); eps = 1 / L
D2 = [H[k + 2] - 2 * H[k + 1] + H[k] for k in range(4)]
fU = lambda U: mp.log1p(1 / U)
def resid(y):
    a, t, w, de = y
    A = n0 / (a * eps); C = t * n0; B = (t + w * eps) * n0; D = n0 + mp.mpf(1) / 2 + de * eps * n0
    F = [fU(A + k) + fU(C + k) - fU(B + k) + fU(n0 + mp.mpf(1) / 2 + k) - fU(D + k) - D2[k] for k in range(4)]
    # rescale: E_j = Delta^j F(0) * n^{j+1} L / ((-1)^j j! (j+1))  (O(1) equations, same zero set)
    E = []
    for j in range(4):
        dj = sum((-1) ** (j - i) * mp.binomial(j, i) * F[i] for i in range(j + 1))
        E.append(dj * n0 ** (j + 1) * L / ((-1) ** j * mp.factorial(j) * (j + 1)))
    return E
y = [mp.mpf(1) / 3, mp.mpf(2), mp.mpf(16) / 3, mp.mpf(1) / 3]
for it in range(40):
    E = resid(y); hstep = mp.mpf(10) ** -80
    J = mp.matrix(4, 4)
    for c in range(4):
        yp = list(y); yp[c] += hstep; Ep = resid(yp)
        for r_ in range(4): J[r_, c] = (Ep[r_] - E[r_]) / hstep
    dy = mp.lu_solve(J, mp.matrix(E))
    y = [y[i] - dy[i] for i in range(4)]
    if max(abs(dy[i]) for i in range(4)) < mp.mpf(10) ** -150: break
a, t, w, de = y
A = n0 / (a * eps); C = t * n0; B = (t + w * eps) * n0; D = n0 + mp.mpf(1) / 2 + de * eps * n0
print("n0 = 1e%d, L = %s, eps = %s, cell index (L*=%d, %d cells) = %d" % (e10, mp.nstr(L, 12), mp.nstr(eps, 12), LS, NE, int(eps * LS * NE)))
print("Newton: %d iterations, residual %s" % (it + 1, mp.nstr(max(abs(v) for v in resid(y)), 3)))
print("exact matching (a, t, w, delta) = (%s, %s, %s, %s), B/n = %s, t - D/n = %s" % tuple(mp.nstr(v, 12) for v in (a, t, w, de, B / n0, t - D / n0)))
X = [n0 ** (j + 1) * L * sum((-1) ** (j + 2 - i) * mp.binomial(j + 2, i) * H[i] for i in range(j + 3)) for j in range(4)]
print("X_j = n^{j+1} L Delta^{j+2} h(n0): %s" % ", ".join(mp.nstr(v, 15) for v in X))
h6 = sum((-1) ** (6 - i) * mp.binomial(6, i) * H[i] for i in range(7))     # = h^(6) at a point of [n0, n0+6]
z = 3
E6 = h6 - mp.psi(5, n0 + mp.mpf(1) / 2 + z) - mp.psi(5, A + z) - mp.psi(5, C + z) + mp.psi(5, B + z) + mp.psi(5, D + z)
print("n0^5 L E_F^(6) (at z ~ 3, via Delta^6 h) = %s" % mp.nstr(n0 ** 5 * L * E6, 12))
print("n0^5 L h^(6) part = %s (leading Q6 = 48 + 26t - 28t^2.., t = 2 eps: %s)" % (mp.nstr(n0 ** 5 * L * h6, 12), mp.nstr(48 + 52 * eps - 112 * eps ** 2, 12)))
print("time %.1f s" % (time.time() - t0))
