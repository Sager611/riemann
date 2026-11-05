# T3: Huisken/Colding-Minicozzi-type monotone functional on the complex side (the brief's bridge):
#   S(t) = sum_{Im z_j > 0} Im z_j,  V(t) = sum_{Im z_j > 0} (Im z_j)^2  -- nonincreasing, zero iff real-rooted.
# (A) its "one-dimensional" representation is Littlewood/Jensen: S = n*Y/2 - (1/2pi) int [log|p(x+iY)| - log|p(x)|] dx (polynomial, Y above all zeros)
# (B) Davenport-Heilbronn test: track the DH off-line zero 0.808517 + 85.699348 i under H_t^DH, compare with a planted zeta pair.
import json, os, time
import numpy as np, mpmath as mp
t00 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__)); SCR = os.path.join(HERE, "..", "..", "..", ".scratch", "explore7")
gam = json.load(open(os.path.join(HERE, "zeros50.json"))); xs = 2 * np.array(gam)   # zeros50.json is written by t1_stam.py
# ---------- (A) Littlewood check on the planted truncated-Xi polynomial ----------
mp.mp.dps = 20
k, Y = 24, 0.30
a, b = xs[k], xs[k + 1]; m = 0.5 * (a + b)
pos = [v for v in xs if v not in (a, b)]
roots = [mp.mpf(v) for v in pos] + [mp.mpf(-v) for v in pos] + [mp.mpc(m, Y), mp.mpc(m, -Y), mp.mpc(-m, Y), mp.mpc(-m, -Y)]
n = len(roots); Ytop = mp.mpf(1)
def integrand(x):
    return sum(mp.log(abs(x + 1j * Ytop - r)) - mp.log(abs(x - r)) for r in roots)
br = sorted(set([float(mp.re(r)) for r in roots]))
pts = [-4000] + br + [4000]
I = mp.quad(integrand, pts)
tail = 2 * sum(0 for r in roots)   # tails: integrand ~ n*Ytop^2/(2x^2) beyond |x|=4000
tail = n * Ytop**2 / 4000          # int_{|x|>L} n Y^2/(2x^2) dx = n Y^2 / L
S_num = n * Ytop / 2 - (I + tail) / (2 * mp.pi)
print(f"(A) planted polynomial n={n}: S true = {2*Y:.6f}; Littlewood S = nY/2 - I/2pi = {mp.nstr(S_num,8)}  (n*Y/2 = {n*float(Ytop)/2}, I/2pi = {mp.nstr((I+tail)/(2*mp.pi),10)})")

# ---------- (B) Davenport-Heilbronn heat flow ----------
mp.mp.dps = 25
kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
c1 = (1 - 1j * kap) / 2; c2 = mp.conj(c1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def xiDH(s):
    return (5 / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2) * (c1 * mp.dirichlet(s, chi) + c2 * mp.dirichlet(s, chib))
def H0DH(z):
    return xiDH(mp.mpf(0.5) + 1j * z / 2) / 8
def xi(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
def H0Z(z):
    return xi(mp.mpf(0.5) + 1j * z / 2) / 8
uG, wG = np.polynomial.hermite.hermgauss(28)
def heat(F, t):
    st = mp.sqrt(t)
    return lambda z: sum(mp.mpf(w) * F(z + 2j * st * mp.mpf(u)) for u, w in zip(uG, wG)) / mp.sqrt(mp.pi) if t > 0 else F(z)
def secant(f, z0, tol=1e-14, it=40):
    za, zb = mp.mpc(z0), mp.mpc(z0) + mp.mpc(1e-5, 1e-5)
    fa, fb = f(za), f(zb)
    for _ in range(it):
        if fb == fa: break
        zc = zb - fb * (zb - za) / (fb - fa)
        za, fa, zb, fb = zb, fb, zc, f(zc)
        if abs(zb - za) < tol: break
    return zb
zDH = mp.mpc(2 * 85.699348485378, 2 * (0.808517182457 - 0.5))   # upper member of the conjugate pair in z = -2i(s-1/2) coordinates (mirror)
fz = lambda z: H0DH(z)
z0 = secant(fz, zDH)
print(f"(B) DH off-line zero at t=0: z = {mp.nstr(z0,14)}  (s = {mp.nstr(0.5 + 1j*mp.conj(z0)/2,14)}), |H0DH(z0)|/|H0DH(z0+1)| = {mp.nstr(abs(fz(z0))/abs(fz(z0+1)),3)}")
# validation of the Gauss-Hermite heat flow: 28 vs 40 nodes at t = 0.1
uG40, wG40 = np.polynomial.hermite.hermgauss(40)
def heat40(F, t):
    st = mp.sqrt(t)
    return lambda z: sum(mp.mpf(w) * F(z + 2j * st * mp.mpf(u)) for u, w in zip(uG40, wG40)) / mp.sqrt(mp.pi)
zt = z0 + mp.mpc(0.3, -0.2)
v28, v40 = heat(H0DH, 0.1)(zt), heat40(H0DH, 0.1)(zt)
print(f"    GH check t=0.1: rel diff 28 vs 40 nodes = {mp.nstr(abs(v28-v40)/abs(v40),3)}")
def track(F, zstart, label, dt=0.01, ystop=0.07):
    z = zstart; t = 0.0; rows = [(0.0, float(mp.im(z)))]
    while True:
        t = round(t + dt, 10)
        z = secant(heat(F, t), z)
        y = float(mp.im(z)); rows.append((t, y))
        if y < ystop or t > 1.0: break
    # extrapolate y^2 -> 0 linearly over the last 3 points before landing (y > 0.03); d(y^2)/dt -> -2 as y -> 0
    pre = [r for r in rows if r[1] > 0.03]
    T = np.array([r[0] for r in pre[-3:]]); Y2 = np.array([r[1] ** 2 for r in pre[-3:]])
    c = np.polyfit(T, Y2, 1); tstar = -c[1] / c[0]
    y0 = rows[0][1]
    print(f"    {label}: y(0)={y0:.6f}, de Bruijn y0^2/2 = {y0*y0/2:.6f}, landing t* = {tstar:.5f} (gain {100*(1-tstar/(y0*y0/2)):.1f}%), x drift {float(mp.re(z)-mp.re(zstart)):+.4f}; y(t): " +
          ", ".join(f"{r[0]:.2f}:{r[1]:.4f}" for r in rows[::3]))
    return tstar, rows
tDH, rowsDH = track(H0DH, z0, "DH actual off-line pair")
# zeta with the same pair planted between the neighbouring zeta zeros (gamma = 84.7355, 87.4253)
gl, gr = 2 * mp.mpf("84.735492980517"), 2 * mp.mpf("87.425274613125")
w = z0
FZ = lambda z: H0Z(z) * (z - w) * (z - mp.conj(w)) / ((z - gl) * (z - gr))
zz0 = secant(FZ, w)
tZ, rowsZ = track(FZ, zz0, "zeta planted pair")
# monotone functional for the DH pair (only one upper zero in the window): S = y, V = y^2 -> check monotone
ys = [r[1] for r in rowsDH]
print(f"    DH: S(t)=y(t) monotone decreasing: {all(ys[i+1] < ys[i] for i in range(len(ys)-1))} -> Lambda_DH >= {tDH:.4f} > 0")
print(f"time {time.time()-t00:.1f}s")
