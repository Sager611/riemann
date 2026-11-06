# Sanity checks for the six-match chain (NOT part of the proof).  uv run --no-project --with mpmath python3 sanity6.py
import mpmath as mp, sys
mp.mp.dps = 300
# ---------- (a) model polynomial at moderate size: localization and critical radius ----------
def pF_coeffs(d, A, B, C, D):
    c = [mp.mpf(1)]
    for k in range(d):
        c.append(c[-1] * D * (k - d) * (k + A) * (k + C) / (A * C * (k + 1) * (k + B) * (k + D)))
    return c
def polyder(c): return [j * c[j] for j in range(1, len(c))]
def peval(c, y): return mp.polyval(c[::-1], y)
print("(a) model p_F: roots within Cloc sqrt(Bd) of B (certified Cloc <= 4.0006) and critical radius max_k |T_k|^(1/k) / sqrt(Bd) (certified <= 2.3)")
for (n, d, eps, (a, t, w, de)) in [(mp.mpf(10) ** 8, 12, mp.mpf(1) / 80, (0.3235, 1.9557, 5.180, 0.3419)),
                                   (mp.mpf(10) ** 9, 20, mp.mpf(1) / 80, (0.3235, 1.9557, 5.180, 0.3419)),
                                   (mp.mpf(10) ** 10, 30, mp.mpf(1) / 200, (0.330, 1.983, 5.27, 0.337))]:
    a, t, w, de = [mp.mpf(x) for x in (a, t, w, de)]
    A = n / (a * eps); C = t * n; B = (t + w * eps) * n; D = n + mp.mpf(1) / 2 + de * eps * n
    c = [cj * B ** j for j, cj in enumerate(pF_coeffs(d, A, B, C, D))]     # variable x = y/B
    sg = mp.sqrt(d / B)                                                   # shift x = 1 + sg*u (well conditioned)
    q = [mp.mpf(0)] * (d + 1)
    for j, cj in enumerate(c):                                            # q(u) = p(1 + sg u)
        for i in range(j + 1):
            q[i] += cj * mp.binomial(j, i) * sg ** i
    roots = mp.polyroots(q[::-1], maxsteps=4000, extraprec=2000)
    imx = max(abs(mp.im(r)) for r in roots)
    assert imx < mp.mpf(10) ** -20, imx
    ru = sorted(mp.re(r) for r in roots)
    loc = max(abs(u) for u in ru) * sg * B / mp.sqrt(B * d)
    dq = polyder(q)
    cp = mp.polyroots(dq[::-1], maxsteps=4000, extraprec=2000)
    ders = [q]
    for k in range(d): ders.append(polyder(ders[-1]))
    worst = 0
    for uu in cp:
        uu = mp.re(uu); q0 = peval(q, uu); x = 1 + sg * uu
        for k in range(1, d + 1):
            Tk = x ** k * sg ** (-k) * peval(ders[k], uu) / q0
            worst = max(worst, abs(Tk) ** (mp.mpf(1) / k) / mp.sqrt(B * d))
    rr = ru
    print("   n=%s d=%d eps=%s: max|root-B|/sqrt(Bd) = %.4f, max_k |T_k|^(1/k)/sqrt(Bd) = %.4f, simple positive roots: %s" % (
        mp.nstr(n, 3), d, mp.nstr(eps, 3), loc, worst, all(rr[i] < rr[i + 1] for i in range(d - 1)) and 1 + sg * rr[0] > 0))
# ---------- (b) residual at the wedge boundary, from the saddle main term ----------
mp.mp.dps = 420
def u0(z):
    u = mp.log(z) / 2 if not isinstance(z, mp.mpc) else mp.log(z) / 2
    u = mp.findroot(lambda v: v * (mp.pi * mp.exp(2 * v) - mp.mpf(9) / 4) - z, mp.re(mp.log(z)) / 2 - mp.log(mp.log(abs(z))) / 2 if True else u)
    return u
def u0n(z, guess):
    return mp.findroot(lambda v: v * (mp.pi * mp.exp(2 * v) - mp.mpf(9) / 4) - z, guess)
def hmain(z, guess):           # Psi(z) - (1/2) log(-f''(u0)) : main part of h = log Mt (G minus log(1+delta))
    u = u0n(z, guess)
    f = 2 * z * mp.log(u) + 9 * u / 2 - mp.pi * mp.exp(2 * u)
    fpp = -2 * z / u ** 2 - 4 * mp.pi * mp.exp(2 * u)
    return f - mp.log(-fpp) / 2, u
n = mp.mpf(10) ** 40
K = 24
d = mp.floor(mp.cbrt(n ** 2 * mp.log(n + 2) / K))
ug = mp.findroot(lambda v: v * (mp.pi * mp.exp(2 * v) - mp.mpf(9) / 4) - n, mp.mpf(44))
L = 2 * ug; eps = 1 / L
ell_cache = {}
def ell(z):                     # log gamma(z) up to const + linear:  h(z) + log Gamma(z+1) - log Gamma(2z+1)
    hz, _ = hmain(z, ug)
    return hz + mp.loggamma(z + 1) - mp.loggamma(2 * z + 1)
lv = [ell(n + k) for k in range(6)]
F = [lv[k + 2] - 2 * lv[k + 1] + lv[k] for k in range(4)]          # Delta^2 log gamma(n+k), k = 0..3
fU = lambda U, k: mp.log(U + k + 1) - mp.log(U + k)
def eqs(a, t, w, de):
    A = n / (a * eps); C = t * n; B = (t + w * eps) * n; D = n + mp.mpf(1) / 2 + de * eps * n
    g = [F[k] - (fU(A, k) + fU(C, k) - fU(B, k) - fU(D, k)) for k in range(4)]
    dg = [g[0], g[1] - g[0], g[2] - 2 * g[1] + g[0], g[3] - 3 * g[2] + 3 * g[1] - g[0]]     # Delta^j at 0
    return [dg[j] * n ** (j + 1) * L for j in range(4)]
sol = mp.findroot(eqs, (mp.mpf('0.33'), mp.mpf('1.98'), mp.mpf('5.27'), mp.mpf('0.337')))
a, t, w, de = sol
A = n / (a * eps); C = t * n; B = (t + w * eps) * n; D = n + mp.mpf(1) / 2 + de * eps * n
print("(b) n = 1e40 (L = %s), d = %s on the wedge boundary n^2 log(n+2) = 24 d^3; branch (a,t,w,delta) = (%s, %s, %s, %s)" % (
    mp.nstr(L, 8), mp.nstr(d, 6), mp.nstr(a, 8), mp.nstr(t, 8), mp.nstr(w, 8), mp.nstr(de, 8)))
def mu(z): return (mp.loggamma(A + z) - mp.loggamma(A) + mp.loggamma(C + z) - mp.loggamma(C) - mp.loggamma(B + z) + mp.loggamma(B)
                   - mp.loggamma(D + z) + mp.loggamma(D))
l0 = lv[0]; s1 = lv[1] - l0 - mu(1)
def EF(z): return ell(n + z) - l0 - mu(z) - z * s1
print("   node check E_F(0..5) =", [mp.nstr(abs(EF(j)), 3) for j in range(6)])
Kr, lam = mp.mpf('2.3'), mp.mpf('1.33')
r = Kr * mp.sqrt(B * d); R = lam * r
best = 0
for x in [d * mp.mpf(i) / 8 for i in range(9)]:
    for th in [mp.pi * j / 12 for j in range(-12, 13)]:
        z = x + R * mp.expj(th)
        v = abs(EF(z)); best = max(best, v)
z = d + R
print("   max |E_F| on the boundary of Omega (sampled) = %s; |E_F(d + lambda r)| = %s" % (mp.nstr(best, 6), mp.nstr(abs(EF(z)), 6)))
bound = mp.mpf('18.752') * d ** 3 / (n ** 2 * L)
print("   certified bound Cfin d^3/(n^2 L) = %s; requirement thetaE = 0.86329" % mp.nstr(bound, 6))
# E_F(6) = Delta^6 E_F(0) = E_F^(6)(xi), xi in (0,6), since E_F(0..5) = 0
print("   n^5 L E_F(6) = n^5 L E_F^(6)(xi) = %s  (certified sup over the disk: <= 2.0504; leading value -2)" % mp.nstr(n ** 5 * L * EF(6), 8))
