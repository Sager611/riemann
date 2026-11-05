# T2: Li-Yau / Hamilton (log-concavity, Harnack) -> Turan/Laguerre.  Questions:
#  (a) is the PF2 / log-concavity input of the Turan inequalities sensitive to t (i.e. to Lambda)?
#  (b) Davenport-Heilbronn: does its kernel have the same positivity / log-concavity?
#  (c) exact Laguerre-hierarchy identity d/dt L_k = (k+1)(2k+1) L_{k+1} - L_k''/2 (sanity check).
import time, mpmath as mp
t0 = time.time()
mp.mp.dps = 30
pi = mp.pi
def Phi(u):          # Riemann kernel, H_0(z) = int_0^inf Phi(u) cos(zu) du = Xi(z/2)/8
    u = mp.mpf(u); s = 0
    for n in range(1, 12):
        e = pi * n * n * mp.e ** (4 * u)
        if e > 400: break
        s += (2 * pi**2 * n**4 * mp.e ** (9 * u) - 3 * pi * n * n * mp.e ** (5 * u)) * mp.e ** (-e)
    return s
# Davenport-Heilbronn: f = (1-ik)/2 L(s,chi) + (1+ik)/2 L(s,chibar), chi mod 5, chi(2) = i
kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
c1 = (1 - 1j * kap) / 2
chi = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
a = lambda n: 2 * mp.re(c1 * chi[n % 5])
def thetaDH(x):
    s = 0
    for n in range(1, 80):
        e = pi * n * n * x / 5
        if e > 300: break
        s += a(n) * n * mp.e ** (-e)
    return s
def PhiDH(v):        # H_0^DH(z) = Xi_DH(z/2)/8 = int_0^inf PhiDH(v) cos(zv) dv
    v = mp.mpf(v)
    return mp.e ** (3 * v) * thetaDH(mp.e ** (4 * v))
print("kappa =", mp.nstr(kap, 15), " a(1..10) =", [mp.nstr(a(n), 6) for n in range(1, 11)])
print("DH kernel evenness Phi(v)/Phi(-v):", [mp.nstr(PhiDH(v) / PhiDH(-v), 12) for v in (0.05, 0.2, 0.5)])
print("zeta kernel evenness:", [mp.nstr(Phi(v) / Phi(-v), 12) for v in (0.05, 0.2, 0.5)])

# (b) positivity and log-concavity of K(sqrt w) on w > 0 (Csordas-Norfolk-Varga input)
def lc2(F, w):  # second derivative of log F(sqrt w)
    return mp.diff(lambda ww: mp.log(F(mp.sqrt(ww))), w, 2)
grid = [mp.mpf(i) / 40 for i in range(1, 161)]   # w in (0, 4]  (u up to 2)
for name, F in (("zeta", Phi), ("DH", PhiDH)):
    vals = [F(mp.sqrt(w)) for w in grid]
    neg = [float(w) for w, v in zip(grid, vals) if v <= 0]
    l2 = [lc2(F, w) for w in grid]
    bad = [float(w) for w, d in zip(grid, l2) if d >= 0]
    print(f"(b) {name}: min kernel on u in (0,2] = {mp.nstr(min(vals),5)} (nonpositive at {neg[:5]}); "
          f"max d2/dw2 log K(sqrt w) = {mp.nstr(max(l2),5)}; log-concavity fails at w = {bad[:6]}")

# (a) Turan inequalities b_m^2 > (2m-1)/(2m+1) b_{m-1} b_{m+1} for H_t, t in a range (Lambda >= 0 => non-real zeros for t < 0)
def moments(F, t, M):
    return [mp.quad(lambda u: u ** (2 * m) * mp.e ** (t * u * u) * F(u), [0, 0.25, 0.5, 1, 1.5, 2.5, 4]) for m in range(M + 2)]
M = 25
for name, F in (("zeta", Phi), ("DH", PhiDH)):
    for t in (-3, -1, -0.3, 0, 0.3):
        b = moments(F, t, M)
        r = [b[m] ** 2 / (b[m - 1] * b[m + 1]) * (2 * m + 1) / (2 * m - 1) for m in range(1, M + 1)]
        print(f"(a) {name} t={t:5}: min_m Turan ratio (>1 means holds) = {mp.nstr(min(r),8)} at m={1+r.index(min(r))}; m=1..4: {[mp.nstr(x,6) for x in r[:4]]}")

# (c) Laguerre hierarchy identity on a small complex-rooted polynomial (exact algebra; numeric check)
mp.mp.dps = 40
roots = [mp.mpf(v) for v in (-7.3, -4.1, -1.2, 2.2, 5.9, 8.8)] + [mp.mpc(0.7, 0.9), mp.mpc(0.7, -0.9)]
coef = [mp.mpf(1)]
for r_ in roots:
    coef = [ (coef[i] if i < len(coef) else 0) - r_ * (coef[i - 1] if i >= 1 else 0) for i in range(len(coef) + 1)]
coef = [mp.re(c) for c in coef]       # p(x) = sum coef[i] x^i, real
def heat(cf, t):   # exp(-t D^2) p
    out = [mp.mpf(0)] * len(cf); d = list(cf); k = 0; fac = mp.mpf(1)
    while any(d):
        for i in range(len(d)): out[i] += (-t) ** k / fac * d[i]
        d = [d[i + 2] * (i + 2) * (i + 1) for i in range(len(d) - 2)] + [0, 0]
        k += 1; fac *= k
    return out
def L(cf, k, x):   # generalized Laguerre: coefficient of y^(2k) in |p(x+iy)|^2
    der = [mp.polyval(list(reversed(cf)), x)]; d = list(cf)
    for j in range(1, 2 * k + 1):
        d = [d[i + 1] * (i + 1) for i in range(len(d) - 1)] or [0]
        der.append(mp.polyval(list(reversed(d)), x))
    return sum((-1) ** (j + k) * der[j] * der[2 * k - j] / (mp.factorial(j) * mp.factorial(2 * k - j)) for j in range(2 * k + 1))
x0, t1, h = mp.mpf("0.6"), mp.mpf("0.05"), mp.mpf("1e-12")
for k in (0, 1, 2):
    lhs = (L(heat(coef, t1 + h), k, x0) - L(heat(coef, t1 - h), k, x0)) / (2 * h)
    Lk = lambda xx: L(heat(coef, t1), k, xx)
    rhs = (k + 1) * (2 * k + 1) * L(heat(coef, t1), k + 1, x0) - mp.diff(Lk, x0, 2) / 2
    print(f"(c) k={k}: dL_k/dt = {mp.nstr(lhs,12)}, (k+1)(2k+1)L_(k+1) - L_k''/2 = {mp.nstr(rhs,12)}")
print(f"time {time.time()-t0:.1f}s")
